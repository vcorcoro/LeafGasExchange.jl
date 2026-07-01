module LeafGasExchange

using Cropbox

include("vaporpressure.jl")
include("weather.jl")
include("diffusion.jl")
include("nitrogen.jl")
include("base.jl")
include("c3.jl")
include("c4.jl")
include("cam.jl")
include("boundarylayer.jl")
include("stomata.jl")
include("intercellularspace.jl")
include("irradiance.jl")
include("energybalance.jl")
include("zeitgeibertime.jl")

@system ModelBase(
    Weather, Nitrogen,
    BoundaryLayer, StomataBase, IntercellularSpace, Irradiance, EnergyBalance,
    ZeitgeiberTime
)

@system ModelBaseDyn(
    Weather, Nitrogen, BoundaryLayer, StomataDyn, IntercellularSpace, 
    Irradiance, EnergyBalanceDyn, ZeitgeiberTime
)

@system ModelC3BB(ModelBase, StomataBallBerry, C3, Controller)
@system ModelC4BB(ModelBase, StomataBallBerry, C4, Controller)

@system ModelC3MD(ModelBase, StomataMedlyn, C3, Controller)
@system ModelC4MD(ModelBase, StomataMedlyn, C4, Controller)

# TODO: include Controller/MinuteController as mixin or not? 
@system ModelCAMBB(ModelBase, StomataBallBerry, CAM)
@system ModelCAMMD(ModelBase, StomataMedlyn, CAM)
@system ModelCAMKB(ModelBase, CAMDyn, StomataKirschbaumCAM, IntercellularSpaceDynCAM)
@system ModelCAMDyn(ModelBaseDyn, CAMDyn, IntercellularSpaceDynCAM)

@system ModelC3KB(ModelBase, StomataKirschbaumC3, IntercellularSpaceDynC3, C3)
@system ModelC3Dyn(ModelBaseDyn, StomataDyn, IntercellularSpaceDynC3, C3)

@system ModelC4Dyn(ModelBaseDyn, StomataDyn, IntercellularSpaceDynC4, C4Dyn)

export ModelC3BB, ModelC3MD, ModelC4BB, ModelC4MD, ModelCAM, ModelCAMKB, ModelC3KB, ModelCAMDyn, ModelC3Dyn, ModelC4Dyn

include("canopy.jl")

include("precompile.jl")

end
