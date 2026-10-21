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
    Weather, Nitrogen, BoundaryLayer, #=StomataDynBase,=# IntercellularSpace, 
    Irradiance, EnergyBalanceTE, ZeitgeiberTime
)

@system ModelC3BB(ModelBase, StomataBallBerry, C3, Controller)
@system ModelC4BB(ModelBase, StomataBallBerry, C4, Controller)

@system ModelC3MD(ModelBase, StomataMedlyn, C3, Controller)
@system ModelC4MD(ModelBase, StomataMedlyn, C4, Controller)

@system ModelCAM0(ModelBase, StomataCAM0, CAM0, Controller)

@system ModelCAMBB(ModelBase, StomataBallBerry, CAM, Controller)
@system ModelCAMMD(ModelBase, StomataMedlyn, CAM, Controller)
# @system ModelCAMKB(ModelBase, CAMDyn, StomataKirschbaumCAM, IntercellularSpaceDynCAM)
@system ModelCAMDynBB(ModelBaseDyn, CAMDyn, IntercellularSpaceDynCAM, StomataDynBB, Controller)
@system ModelCAMDynMD(ModelBaseDyn, CAMDyn, IntercellularSpaceDynCAM, StomataDynMD, Controller)

# @system ModelC3KB(ModelBase, StomataKirschbaumC3, IntercellularSpaceDynC3, C3)
@system ModelC3DynBB(ModelBaseDyn, StomataDynBB, IntercellularSpaceDynC3, C3, Controller)
@system ModelC3DynMD(ModelBaseDyn, StomataDynMD, IntercellularSpaceDynC3, C3, Controller)

@system ModelC4DynBB(ModelBaseDyn, StomataDynBB, IntercellularSpaceDynC4, C4Dyn, Controller)
@system ModelC4DynMD(ModelBaseDyn, StomataDynMD, IntercellularSpaceDynC4, C4Dyn, Controller)

@system ModelC4ExactStsMD(ModelBase, StomataMedlyn, C4Exact, Controller)
@system ModelC4ExactStsBB(ModelBase, StomataBallBerry, C4Exact, Controller)

@system ModelC4ExactDynBB(ModelBaseDyn, StomataDynBB, IntercellularSpace, C4Exact, Controller)
@system ModelC4ExactDynMD(ModelBaseDyn, StomataDynMD, IntercellularSpace, C4Exact, Controller)

@system ModelC3StsMD(ModelBase, StomataMedlyn, C3, Controller)
@system ModelC3StsBB(ModelBase, StomataBallBerry, C3, Controller)
@system ModelC3StsMDTE(ModelBase, StomataMedlyn, EnergyBalanceTE, C3, Controller)
@system ModelC3StsBBTE(ModelBase, StomataBallBerry, EnergyBalanceTE, C3, Controller)
@system ModelC3DynMD(ModelBaseDyn, StomataDynMD, IntercellularSpaceDynC3, C3, Controller)
@system ModelC3DynBB(ModelBaseDyn, StomataDynBB, IntercellularSpaceDynC3, C3, Controller)

# export ModelC3BB, ModelC3MD, ModelC4BB, ModelC4MD, ModelCAMBB, ModelC3KB, ModelCAMDynBB, ModelCAMDynMD, ModelC3Dyn, ModelC4Dyn
# export ModelC3BB, ModelC3MD, ModelC3DynBB, ModelC3DynMD
export ModelC3StsMD, ModelC3StsBB, ModelC3StsMDTE, ModelC3StsBBTE, ModelC3DynMD, ModelC3DynBB
export ModelC4BB, ModelC4MD, ModelC4DynBB, ModelC4DynMD
export ModelCAMBB, ModelCAMMD, ModelCAMDynBB, ModelCAMDynMD, ModelCAM0
export ModelC4ExactStsBB, ModelC4ExactStsMD, ModelC4ExactDynBB, ModelC4ExactDynMD




include("canopy.jl")

include("precompile.jl")

end
