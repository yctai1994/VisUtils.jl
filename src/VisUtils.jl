module VisUtils

import FFTW
import DelimitedFiles as dlm

const VecIO = AbstractVector
const MatIO = AbstractMatrix

include("./StatsVector.jl")
include("./TimeDelay.jl")
include("./Spectrum.jl")

end # module VisUtils
