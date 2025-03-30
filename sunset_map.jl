#!/bin/julia

import ThreadsX

using Dates
using Insolation
import Insolation.Parameters as IP
import Insolation.OrbitalData
import ClimaParams as CP

# adapted from Insolation tests
param_set = IP.InsolationParameters(Float32) # ClimaParams has a side effect here, that's fun
rtol = 1e-2
od = Insolation.OrbitalData()

date0 = DateTime("2000-01-01T11:58:56.816")

function below_horizon(date, lon, lat)
    args = (Insolation.helper_instantaneous_zenith_angle(date, date0, od, param_set)...,lon, lat,)
    sza, azi, d = instantaneous_zenith_angle(args...)
    return (pi/2 - sza, azi) # height above horizon and e/w direction: azi approx = 2pi => east, approx 1pi => west
end

# cool so for sunset, sza ~ approx pi/2 and azi ~ approx pi
function sunset(pos; date=DateTime(2020, 2, 1, 0))
    argmin(t -> begin
        horiz, azi = below_horizon(t, pos[1], pos[2]) # lon, lat
        return abs(horiz) + (azi > 1.5 * pi ? 1000.0 : 0.0) # avoid sunrise, lol
    end, date .+ (Minute(0):Minute(1):Minute(60*24)))
end

# lon/lat grid in degrees
lons = Float32.(-180:5:180)
lats = Float32.(-90:5:90)
positions = Iterators.product(lons, lats)

sunsets = ThreadsX.map(sunset, positions) # todo - work out what to do when there is no sunset (currently we're just using lowest point in the west. which is probably fine)

# plotting
using GLMakie, GeoMakie # cairomakie didn't work for me :(

GLMakie.activate!()

fig = Figure()
ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    dest = "+proj=wintri", # the CRS in which you want to plot
)
# You can plot your data the same way you would in Makie
# scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2))
# fig

# need a circular colour map where 0:1 are the same
cmap = to_colormap(:plasma)

cmap_circ = [cmap[1:2:256]...; cmap[256:-2:1]...]

# s = surface!(ga, lons, lats, map(x -> Time(x).instant.value, sunsets) ./ Time(23,59,59,999).instant.value;  colormap = cmap_circ, shading = NoShading)
nsfrac2time = x -> string(Time(Nanosecond(round(x * Time(23,59,59,999).instant.value))))
Colorbar(fig[1,2], s; tickformat = v -> nsfrac2time.(v))

lines!(ga, GeoMakie.coastlines(), overdraw=true, color=:black) # plot coastlines from Natural Earth as a reference

# animation
dates = DateTime(2020, 1, 1, 0):Day(1):DateTime(2021, 1, 1, 0)
positionstimes = Iterators.product(lons, lats, dates)

sunsets = ThreadsX.map(x->sunset([x[1], x[2]], date=x[3]), positionstimes) # todo - work out what to do when there is no sunset (currently we're just using lowest point in the west. which is probably fine)

frames = size(sunsets)[3]
# record(fig, "sunset_map.mp4", 1:frames; framerate=10) do i
# end
#for i in 1:frames
record(fig, "sunset_map.mp4", 1:frames; framerate=24) do i
    surface!(ga, lons, lats, map(x -> Time(x).instant.value, sunsets[:, :, i]) ./ Time(23,59,59,999).instant.value;  colormap = cmap_circ, shading = NoShading, overdraw=true)
    lines!(ga, GeoMakie.coastlines(), overdraw=true, color=:black) # plot coastlines from Natural Earth as a reference
    ga.title = string(Date(dates[i]))
    # sleep(0.1)
end # overdraw is important - stops weird z-fighting



# zoomed in on europe
fig = Figure()
ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    dest = "+proj=wintri", # the CRS in which you want to plot
)
i = 355
europe = sunsets[-10 .< lons .< 30, 35 .< lats .< 70, i]
normify = x -> min(max(0, (Time(x) - Time(minimum(europe))).value ./ (Time(maximum(europe)) - Time(minimum(europe))).value), 1)
# max 0,1 is important otherwise makie helpfully descales our clamping

denormify = x -> string(Time(Nanosecond(Time(minimum(europe)).instant.value + (Time(maximum(europe)) - Time(minimum(europe))).value * x)))
ga.title = string(Date(dates[i]))
ylims!(ga, 35, 70)
xlims!(ga, -10, 30)
ilat_min = findfirst(x -> x >= 35, lats)
ilat_max = findlast(x -> x <= 70, lats)
ilon_min = findfirst(x -> x >= -10, lons)
ilon_max = findlast(x -> x <= 30, lons)

s = contour!(ga, lons[ilon_min:ilon_max], lats[ilat_min:ilat_max], map(normify, sunsets[ilon_min:ilon_max, ilat_min:ilat_max, i]);  levels=10, colormap = cmap, overdraw=true, labels=true, labelformatter=v -> denormify.(v))
# s = surface!(ga, lons, lats, map(normify, sunsets[:, :, i]);  colormap = cmap, shading = NoShading, overdraw=true)
lines!(ga, GeoMakie.coastlines(), overdraw=true, color=:black) # plot coastlines from Natural Earth as a reference
# Colorbar(fig[1,2], s; tickformat = v -> denormify.(v), ticks=0:0.1:1)
