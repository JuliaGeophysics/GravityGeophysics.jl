# Author: @pankajkmishra
"""
    Forward.jl - Unified Forward Modeling Interface

Provides a common interface for both Prism and FD forward modeling methods.
Allows method selection via dispatch on ForwardMethod types.
"""

"""
    ForwardMethod

Abstract type for forward modeling methods.
"""
abstract type ForwardMethod end

"""
    Prism <: ForwardMethod

Analytical prism method using closed-form solution.
Exact but O(N_cells × N_receivers) complexity.
"""
struct Prism <: ForwardMethod end

"""
    FiniteDifference <: ForwardMethod

Finite difference method solving Poisson equation.
Scales better for large meshes.

# Fields
- `options::FDOptions`: Solver options
"""
struct FiniteDifference <: ForwardMethod
    options::FDOptions
end

# Default constructor
FiniteDifference() = FiniteDifference(FDOptions())

"""
    forward(model::GravityModel, receivers::GravityReceivers, method::ForwardMethod)

Compute gravity forward response.

# Arguments
- `model`: GravityModel with density distribution
- `receivers`: Observation locations  
- `method`: ForwardMethod (Prism() or FiniteDifference())

# Returns
- `GravityData`: Computed gravity data

# Examples
```julia
# Using Prism method
data = forward(model, receivers, Prism())

# Using FD method with custom options
opts = FDOptions(tol=1e-12, maxiter=8000, use_prism_bc=true)
data = forward(model, receivers, FiniteDifference(opts))
```
"""
function forward(model::GravityModel, receivers::GravityReceivers, method::Prism)
    return forward_prism(model, receivers)
end

function forward(model::GravityModel, receivers::GravityReceivers, method::FiniteDifference)
    data, _ = forward_fd(model, receivers; options=method.options)
    return data
end

# Default to Prism method if no method specified
function forward(model::GravityModel, receivers::GravityReceivers)
    return forward(model, receivers, Prism())
end

"""
    forward_with_info(model::GravityModel, receivers::GravityReceivers, method::ForwardMethod)

Compute gravity forward response with additional solver info.

# Returns
- `GravityData`: Computed gravity data
- `info::Dict`: Solver information
"""
function forward_with_info(model::GravityModel, receivers::GravityReceivers, method::Prism)
    data = forward_prism(model, receivers)
    info = Dict("method" => "Prism", "ncells" => count(!iszero, model.density))
    return data, info
end

function forward_with_info(model::GravityModel, receivers::GravityReceivers, method::FiniteDifference)
    return forward_fd(model, receivers; options=method.options)
end

"""
    compute_stats(data1::GravityData, data2::GravityData)

Compute comparison statistics between two gravity datasets.

# Returns
- `stats::Dict`: Dictionary with RMSE, MAE, MAX error, relative errors
"""
function compute_stats(data1::GravityData, data2::GravityData)
    gz1 = data1.gz
    gz2 = data2.gz
    
    @assert length(gz1) == length(gz2) "Data arrays must have same length"
    
    n = length(gz1)
    s2 = 0.0  # Sum of squared differences
    s1 = 0.0  # Sum of absolute differences
    mx = 0.0  # Max absolute difference
    
    for i in 1:n
        d = gz1[i] - gz2[i]
        s2 += d * d
        s1 += abs(d)
        mx = max(mx, abs(d))
    end
    
    rmse = sqrt(s2 / n)
    mae = s1 / n
    
    # Relative errors
    max_signal = maximum(abs.(gz1))
    rel_rmse = max_signal > 0 ? rmse / max_signal * 100 : 0.0
    rel_max = max_signal > 0 ? mx / max_signal * 100 : 0.0
    
    return Dict(
        "rmse_si" => rmse,
        "mae_si" => mae,
        "max_si" => mx,
        "rmse_mgal" => rmse * SI_TO_MGAL,
        "mae_mgal" => mae * SI_TO_MGAL,
        "max_mgal" => mx * SI_TO_MGAL,
        "rel_rmse_pct" => rel_rmse,
        "rel_max_pct" => rel_max
    )
end

"""
    compare_methods(model::GravityModel, receivers::GravityReceivers;
                    fd_options::FDOptions=FDOptions())

Compare Prism and FD methods on the same problem.

# Returns
- `data_prism::GravityData`: Prism method result
- `data_fd::GravityData`: FD method result
- `stats::Dict`: Comparison statistics
"""
function compare_methods(model::GravityModel, receivers::GravityReceivers;
                         fd_options::FDOptions=FDOptions())
    t0 = time()
    data_prism = forward(model, receivers, Prism())
    t_prism = time() - t0
    
    t0 = time()
    data_fd, fd_info = forward_fd(model, receivers; options=fd_options)
    t_fd = time() - t0
    
    stats = compute_stats(data_prism, data_fd)
    stats["time_prism"] = t_prism
    stats["time_fd"] = t_fd
    stats["fd_iterations"] = fd_info["iterations"]
    
    # Print comparison
    println("═"^60)
    println("Method Comparison Results")
    println("═"^60)
    @printf("  Prism time:     %.3f s\n", t_prism)
    @printf("  FD time:        %.3f s (CG iters: %d)\n", t_fd, fd_info["iterations"])
    println("─"^60)
    @printf("  RMSE:           %.6e m/s² (%.6f mGal)\n", stats["rmse_si"], stats["rmse_mgal"])
    @printf("  MAE:            %.6e m/s² (%.6f mGal)\n", stats["mae_si"], stats["mae_mgal"])
    @printf("  MAX error:      %.6e m/s² (%.6f mGal)\n", stats["max_si"], stats["max_mgal"])
    @printf("  Relative RMSE:  %.3f%%\n", stats["rel_rmse_pct"])
    @printf("  Relative MAX:   %.3f%%\n", stats["rel_max_pct"])
    println("═"^60)
    
    return data_prism, data_fd, stats
end
