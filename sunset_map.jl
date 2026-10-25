#!/usr/bin/env julia

module SunsetMap

using Dates, H3.API, TimeZones
import Arrow, ClimaParams, Insolation, LibGEOS, TimeZoneFinder

const PARAMETERS = Insolation.Parameters.InsolationParameters(Float64)
const ORBIT = Insolation.compute_orbital_parameters(PARAMETERS)
const EPOCH = DateTime(2000, 1, 1, 11, 58, 56, 816)
const HORIZON = sind(-0.833)
# UTC instants from https://www.weather.gov/dvn/Climate_Astronomical_Seasons
const SOLSTICES = (DateTime(2026, 6, 21, 8, 24), DateTime(2026, 12, 21, 20, 50))

function solar_event_utc(date::Date, latitude, longitude, sunrise::Bool)
    midnight = DateTime(date)
    estimate = midnight + Millisecond(round(Int, (720 - 4longitude) * 60_000))
    sinlat, coslat = sincos(deg2rad(latitude))
    direction = sunrise ? -1 : 1
    for _ in 1:12
        _, declination, hourangle = Insolation.distance_declination_hourangle(
            estimate, EPOCH, ORBIT, PARAMETERS, true,
        )
        cosine = (HORIZON - sinlat * sin(declination)) / (coslat * cos(declination))
        abs(cosine) < 1 || return nothing
        uncorrected = 2pi * mod(datetime2julian(estimate), 1)
        equation = mod(hourangle - uncorrected + pi, 2pi) - pi
        minutes = 720 - 4longitude + (direction * acos(cosine) - equation) * 720 / pi
        next = midnight + Millisecond(round(Int, minutes * 60_000))
        abs(Dates.value(next - estimate)) <= 1 && return next
        estimate = next
    end
    name = sunrise ? "Sunrise" : "Sunset"
    error("$name solver failed at ($latitude, $longitude) on $date")
end

sunset_utc(date::Date, latitude, longitude) = solar_event_utc(date, latitude, longitude, false)
sunrise_utc(date::Date, latitude, longitude) = solar_event_utc(date, latitude, longitude, true)

function local_solar_event(date::Date, latitude, longitude, timezone, sunrise::Bool)
    calculate = sunrise ? sunrise_utc : sunset_utc
    event = calculate(date, latitude, longitude)
    isnothing(event) && return nothing
    local_event = astimezone(ZonedDateTime(event, tz"UTC"), timezone)
    # Match the local civil date across the date line.
    if Date(local_event) != date
        event = calculate(date + (date - Date(local_event)), latitude, longitude)
        isnothing(event) && return nothing
        local_event = astimezone(ZonedDateTime(event, tz"UTC"), timezone)
    end
    Date(local_event) == date || return nothing
    return Dates.value(Time(local_event)) / 60_000_000_000
end

local_sunset(date::Date, latitude, longitude, timezone) =
    local_solar_event(date, latitude, longitude, timezone, false)
local_sunrise(date::Date, latitude, longitude, timezone) =
    local_solar_event(date, latitude, longitude, timezone, true)

winter_solstice(latitude) = SOLSTICES[latitude >= 0 ? 2 : 1]

function winter_solar_event(latitude, longitude, timezone, sunrise::Bool)
    solstice = winter_solstice(latitude)
    date = Date(astimezone(ZonedDateTime(solstice, tz"UTC"), timezone))
    calculate = sunrise ? local_sunrise : local_sunset
    return calculate(date, latitude, longitude, timezone)
end

winter_sunset(latitude, longitude, timezone) =
    winter_solar_event(latitude, longitude, timezone, false)
winter_sunrise(latitude, longitude, timezone) =
    winter_solar_event(latitude, longitude, timezone, true)

function global_cells(resolution::Int)
    0 <= resolution <= 15 || throw(ArgumentError("H3 resolution must be in 0:15"))
    cells = UInt64[]
    sizehint!(cells, getNumCells(resolution))
    for base in getRes0Cells()
        append!(cells, cellToChildren(base, resolution))
    end
    return cells
end

split_index(cell::UInt64) = (index_lower=cell % UInt32, index_upper=(cell >> 32) % UInt32)

function timezone_lookup()
    version = TimeZoneFinder._timezone_boundary_builder_version()
    data = TimeZoneFinder.load_data(version)
    context = LibGEOS.GEOSContext()
    meshes = TimeZoneFinder.Meshes
    polygons = map(data.shapes) do shape
        rings = map(meshes.rings(shape.polyarea)) do ring
            points = [collect(meshes.coordinates(p)) for p in meshes.vertices(ring)]
            push!(points, first(points))
        end
        prepared = LibGEOS.prepareGeom(LibGEOS.Polygon(rings, context))
        finalizer(LibGEOS.destroyGeom, prepared)
        prepared
    end
    return (tree=LibGEOS.STRtree(polygons; context), zones=IdDict(zip(polygons, data.tzs)), context)
end

function lookup_timezone(lookup, latitude, longitude)
    point = LibGEOS.Point(longitude, latitude, lookup.context)
    try
        matches = filter(p -> LibGEOS.covers(p, point), LibGEOS.query(lookup.tree, point))
        isempty(matches) && error("Timezone missing at ($latitude, $longitude)")
        # Use a stable choice for overlapping timezone boundaries.
        return lookup.zones[first(sort!(matches; by=p -> string(lookup.zones[p])))]
    finally
        finalize(point)
    end
end

function solar_event_table(cells; event=:sunset)
    event in (:sunrise, :sunset) || throw(ArgumentError("event must be :sunrise or :sunset"))
    calculate = event === :sunrise ? winter_sunrise : winter_sunset
    lookup = timezone_lookup()
    centres = [cellToLatLng(cell)::LatLng for cell in cells]
    # GEOS lookups share one context and run serially. Solar calculations use threads.
    zones = [lookup_timezone(lookup, rad2deg(c.lat), rad2deg(c.lng)) for c in centres]
    values = fill(Float32(Inf), length(cells))
    Threads.@threads for i in eachindex(cells)
        centre = centres[i]
        minutes = calculate(rad2deg(centre.lat), rad2deg(centre.lng), zones[i])
        isnothing(minutes) || (values[i] = minutes)
    end
    keep = findall(isfinite, values)
    return (
        index_lower = [split_index(cells[i]).index_lower for i in keep],
        index_upper = [split_index(cells[i]).index_upper for i in keep],
        value = values[keep],
    )
end

sunset_table(cells) = solar_event_table(cells; event=:sunset)
sunrise_table(cells) = solar_event_table(cells; event=:sunrise)

function generate(output=nothing; resolution=5, event=:sunset)
    event in (:sunrise, :sunset) || throw(ArgumentError("event must be :sunrise or :sunset"))
    output = isnothing(output) ? "$(event)_2026_r$(resolution).arrow" : output
    @info "Enumerating H3 cells" resolution
    cells = global_cells(resolution)
    @info "Calculating local $event times" cells=length(cells) threads=Threads.nthreads()
    table = solar_event_table(cells; event)
    Arrow.write(output, table; compress=nothing, dictencode=false, metadata=Dict(
        "solar_event" => string(event),
        "value_units" => "minutes since local midnight",
        "year" => "2026",
        "h3_resolution" => string(resolution),
        "solstices_utc" => join(string.(SOLSTICES), ","),
        "date_rule" => "local civil date containing each solstice instant",
        "selection" => "December solstice at latitude >= 0; June solstice at latitude < 0",
        "solar_altitude_degrees" => "-0.833",
        "timezone_boundaries" => TimeZoneFinder._timezone_boundary_builder_version(),
        "boundary_source" => "timezone-boundary-builder; OpenStreetMap contributors",
        "boundary_license" => "https://opendatacommons.org/licenses/odbl/1-0/",
        "tzdata" => TimeZones.TZData.tzdata_version(),
        "timezone_overlap_rule" => "lexicographically first timezone identifier",
    ))
    @info "Wrote Arrow file" output rows=length(table.value) omitted=length(cells)-length(table.value)
    return table
end

function main(args)
    sunrise = "--sunrise" in args
    positional = filter(!=("--sunrise"), args)
    length(positional) <= 2 || throw(ArgumentError(
        "Usage: sunset_map.jl [--sunrise] [output.arrow] [resolution]",
    ))
    output = isempty(positional) ? nothing : positional[1]
    resolution = length(positional) == 2 ? parse(Int, positional[2]) : 5
    @time generate(output; resolution, event=sunrise ? :sunrise : :sunset)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    SunsetMap.main(ARGS)
end
