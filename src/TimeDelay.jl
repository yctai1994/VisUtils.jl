export FemtosecondDelay, PicosecondDelay, MillimeterDelay

struct TimeDelayUnit{T} end

const FemtosecondDelay = TimeDelayUnit{:fs}()
const PicosecondDelay  = TimeDelayUnit{:ps}()
const MillimeterDelay  = TimeDelayUnit{:mm}()

export TimeDelay

struct TimeDelay
	vals::Vector{Float64}
	unit::Base.RefValue{TimeDelayUnit}

	function TimeDelay(vals::VecIO{Float64}, unit::TimeDelayUnit)
		return new(vals, Ref{TimeDelayUnit}(unit))
	end
end

export mm2fs, mm2fs!

@inline mm2fs(mm::Real) = mm / 0.000149896229

function mm2fs!(delay::TimeDelay)
    if delay.unit[] ≡ MillimeterDelay
        mm2fs!(delay.vals)
        delay.unit[] = FemtosecondDelay
    end
        return nothing
end

function mm2fs!(delay::VecIO{Float64})
    @simd for i in eachindex(delay)
        @inbounds delay[i] = mm2fs(delay[i])
    end
    return nothing
end

export mm2ps, mm2ps!

@inline mm2ps(mm::Real) = mm / 0.149896229

function mm2ps!(delay::TimeDelay)
    if delay.unit[] ≡ MillimeterDelay
        mm2ps!(delay.vals)
        delay.unit[] = PicosecondDelay
    end
        return nothing
end

function mm2ps!(delay::VecIO{Float64})
    @simd for i in eachindex(delay)
        @inbounds delay[i] = mm2ps(delay[i])
    end
    return nothing
end
