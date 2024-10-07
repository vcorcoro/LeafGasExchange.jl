include("./LeafGasExchange3.jl")
using ..LeafGasExchange3
using Cropbox

# ModelCAM uses StomataCAM and CAM photosynthesis
# otherwise uses same ModelBase as C3 and C4 
Cropbox.hierarchy(LeafGasExchange3.ModelCAM; skipcontext = true)

# system for running with weather data
using CSV
using DataFrames
@system S(LeafGasExchange3.ModelCAM, Controller) begin
    weather_data => CSV.read("sample_data/TempleApril2015Interp30.csv", DataFrame) ~ provide
    T_air(weather_data): air_temperature => weather_data[!, :Temperature] ~ drive(u"°C")
    RH(weather_data): relative_humidity => weather_data[!, :Relative_Humidity] ~ drive(u"percent")
    PHI(weather_data) => weather_data[!, :GHI] ~ drive(u"W/m^2")
    PFD(PHI,Q): photon_flux_density => PHI*Q ~ track(u"μmol/m^2/s")
end

c_weather = @config(
    :Weather => (
        # PFD = 0u"μmol/m^2/s",
        CO2 = 400u"μmol/mol",
        # RH = 60u"percent",
        # T_air = 25u"°C",
        wind = 2.0u"m/s",
    )
)

c = @config(
    c_weather,
    :Clock => :step => 0.5u"hr",
    :TemperatureDependence => :Tb => 20u"°C"
)

s = simulate(S; config=c, stop=7u"d")

# net assimilation
plot(s, :time, [:Asv,:Asc,:A_net]; kind=:line)
# co2_demand_function_ci
plot(s, :time, [:Ac_ci,:Aj_ci,:Ad_ci]; kind=:line)
# co2_demand_function_cc
plot(s, :time, [:Ac_cc,:Aj_cc,:Ad_cc]; kind=:line)

plot(s, :time, :Cc; kind=:line) # corrected_mesophyll_co2
plot(s, :time, :Ci; kind=:line, ylim=(0,450)) # intercellular_co2

plot(s, :time, :gs; kind=:line) # stomatal_conductance_to_water

plot(s, :time, :Tk; kind=:line, ylim=(284,306)) # temperature

plot(s, :time, :E) # transpiration

# circadian oscillator variables: 
plot(s, :time, :m; kind=:line) # malic_acid_concentration
plot(s, :time, :z; kind=:line, ylim=(0,1.2)) # circadian_rhythm_order

# sensitivity analysis of parameters
let y = :A_net,
    group = :0 => :m_max => [2e6, 1.2e8, 2.3e8, 5e8],
    kind = :line
    p = visualize(S, :time, y; group, kind, config=c, stop=1u"d")
    visualize!(p, S, :time, y; kind, config=c, stop=1u"d")
end
