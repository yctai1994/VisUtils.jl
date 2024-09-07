export FrequencyUnit, TeraHertz, PetaHertz

struct FrequencyUnit{T} end

const TeraHertz = FrequencyUnit{:THz}()
const PetaHertz = FrequencyUnit{:PHz}()

export Frequency

struct Frequency
	vals::Vector{Float64}
	unit::Base.RefValue{FrequencyUnit}

	function Frequency(delay::TimeDelay, Δt::Real)
        unit = delay.unit[] ≡ FemtosecondDelay ? PetaHertz :
               delay.unit[] ≡ PicosecondDelay  ? TeraHertz : error("")
        freq = FFTW.rfftfreq(length(delay.vals), inv(Δt))

        return new(freq, Ref{FrequencyUnit}(unit))
    end
end

export real_spectra

function real_spectra(signal::MatIO, background::VecIO; dims::Int = 1)
    nrow, ncol = size(signal)
    nrow ≡ length(background) || error("")

    sub_signal = Matrix{Float64}(undef, nrow, ncol)

    for j in 1:ncol
        @simd for i in eachindex(background)
            @inbounds sub_signal[i,j] = signal[i,j] - background[i]
        end
    end

    return real_spectra(sub_signal; dims = dims)
end

function real_spectra(signal::MatIO; dims::Int = 1)
    rfft = FFTW.rfft(signal, dims)
    spec = Matrix{Float64}(undef, size(rfft)...)
    fact = inv(0.5 * size(signal, dims))

    @inbounds for i in eachindex(spec)
        spec[i] = fact * abs(rfft[i])
    end

    return spec
end
