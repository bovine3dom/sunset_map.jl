using Test
using Dates
using TimeZones
using Arrow
using H3.API
using TimeZoneFinder: timezones_at

include(joinpath(@__DIR__, "..", "sunset_map.jl"))
using .SunsetMap

zone(name) = TimeZone(name, TimeZones.Class(:ALL))
cell(lat, lon) = latLngToCell(LatLng(deg2rad(lat), deg2rad(lon)), 5)

@testset "Solar calculation" begin
    # Independent sea-level reference values from Astral 3.2, in local minutes.
    # https://astral.readthedocs.io/en/latest/
    for (date, lat, lon, timezone, expected) in [
        (Date(2026,12,21), 51.5072, -0.1276, "Europe/London", 953.061),
        (Date(2026,6,21), -33.8688, 151.2093, "Australia/Sydney", 1013.594),
        (Date(2026,12,21), 21.3099, -157.8581, "Pacific/Honolulu", 1074.636),
        (Date(2026,6,21), 61.2181, -149.9003, "America/Anchorage", 1422.035),
        (Date(2026,12,21), 64.1466, -21.9426, "Atlantic/Reykjavik", 928.551),
        (Date(2026,6,21), 1.8721, -157.4278, "Pacific/Kiritimati", 1118.159),
        (Date(2026,6,21), 64.1466, -21.9426, "Atlantic/Reykjavik", 2.650),
    ]
        actual = SunsetMap.local_sunset(date, lat, lon, zone(timezone))
        @test isapprox(actual, expected; atol=1.5)
    end
    @test isnothing(SunsetMap.local_sunset(Date(2026,12,21), 69.6492, 18.9553, zone("Europe/Oslo")))
    @test isnothing(SunsetMap.local_sunset(Date(2026,6,21), 69.6492, 18.9553, zone("Europe/Oslo")))
    @test isnothing(SunsetMap.winter_sunset(-75.0, 0.0, tz"UTC"))
    @test SunsetMap.winter_sunset(64.1466, -21.9426, zone("Atlantic/Reykjavik")) > 900
    @test SunsetMap.winter_sunset(-33.8688, 151.2093, zone("Australia/Sydney")) < 1020
    @test SunsetMap.winter_solstice(0.0) == SunsetMap.SOLSTICES[2]

    # Solstice dates follow civil time on both sides of the date line.
    june, december = ZonedDateTime.(SunsetMap.SOLSTICES, Ref(tz"UTC"))
    @test Date(astimezone(june, zone("Pacific/Honolulu"))) == Date(2026,6,20)
    @test Date(astimezone(december, zone("Pacific/Kiritimati"))) == Date(2026,12,22)
    @test SunsetMap.winter_sunset(1.8721, -157.4278, zone("Pacific/Kiritimati")) ≈
        SunsetMap.local_sunset(Date(2026,12,22), 1.8721, -157.4278, zone("Pacific/Kiritimati"))
end

@testset "H3 indices" begin
    for resolution in (0, 1, 5)
        cells = SunsetMap.global_cells(resolution)
        @test length(cells) == getNumCells(resolution)
        @test allunique(cells)
        @test all(c -> isValidCell(c) && getResolution(c) == resolution, cells)
    end
    @test_throws ArgumentError SunsetMap.global_cells(-1)
    @test_throws ArgumentError SunsetMap.global_cells(16)
    for index in (cell(51.5,-0.1), cell(-33.9,151.2), typemax(UInt64))
        lower, upper = SunsetMap.split_index(index)
        @test lower isa UInt32
        @test upper isa UInt32
        @test UInt64(lower) | (UInt64(upper) << 32) == index
    end
end

@testset "Indexed timezone lookup" begin
    lookup = SunsetMap.timezone_lookup()
    @test string(SunsetMap.lookup_timezone(lookup, -45.5752, -72.0662)) == "America/Coyhaique"
    centres = cellToLatLng.(SunsetMap.global_cells(2))
    for centre in centres
        lat, lon = rad2deg(centre.lat), rad2deg(centre.lng)
        expected = first(sort!(timezones_at(lat, lon); by=string))
        @test SunsetMap.lookup_timezone(lookup, lat, lon) == expected
    end
    for (lat, lon) in [(69.8,-141.0), (50.5,1.0), (52.5061,13.358)]
        expected = first(sort!(timezones_at(lat, lon); by=string))
        @test SunsetMap.lookup_timezone(lookup, lat, lon) == expected
    end
end

@testset "Arrow output" begin
    cells = [cell(51.5072,-0.1276), cell(-33.8688,151.2093), cell(69.6492,18.9553), cell(-75,0)]
    table = SunsetMap.sunset_table(cells)
    @test length(table.value) == 2
    @test all(v -> 0 <= v < 1440, table.value)
    indices = UInt64.(table.index_lower) .| (UInt64.(table.index_upper) .<< 32)
    @test indices == cells[1:2]
    @test isempty(SunsetMap.sunset_table(UInt64[]).value)
    @test isempty(SunsetMap.sunset_table(cells[3:4]).value)
    mktempdir() do dir
        path = joinpath(dir, "sunset.arrow")
        generated = SunsetMap.generate(path; resolution=0)
        restored = Arrow.Table(path)
        @test eltype(restored.index_lower) == UInt32
        @test eltype(restored.index_upper) == UInt32
        @test eltype(restored.value) == Float32
        @test restored.value == generated.value
        @test restored.index_lower == generated.index_lower
        @test restored.index_upper == generated.index_upper
        @test Arrow.getmetadata(restored)["value_units"] == "minutes since local midnight"
        @test all(v -> isfinite(v) && 0 <= v < 1440, restored.value)
    end
end
