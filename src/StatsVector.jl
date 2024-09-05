export StatsVector

# Statistics (vectors) for band plotting
struct StatsVector
	avg::Vector{Float64}
	μpσ::Vector{Float64}
	μmσ::Vector{Float64}

	function StatsVector(dat::MatIO{Float64})
		nrow, ncol = size(dat)

		avg = Vector{Float64}(undef, nrow)
		μpσ = Vector{Float64}(undef, nrow) # μ - σ
		μmσ = Vector{Float64}(undef, nrow) # μ + σ

		# use μpσ as std buffer
		for j in eachindex(1:ncol)
			@inbounds for i in eachindex(avg)
				avg[i], μpσ[i] = welfordStep(avg[i], μpσ[i], dat[i,j], j)
			end
		end

		@simd for i in eachindex(μpσ)
			@inbounds μpσ[i] = sqrt(μpσ[i])
		end

		@inbounds for i in eachindex(avg)
			μmσ[i] = max(0.0, avg[i] - μpσ[i])
			μpσ[i] = μpσ[i] + avg[i]
		end

		return new(avg, μpσ, μmσ)
	end
end
