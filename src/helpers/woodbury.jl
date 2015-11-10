type woodbury_identity
    U::AbstractMatrix
    V::AbstractMatrix
    invA_U::AbstractMatrix
    factored_capacitance_matrix
    factored_A::abstract_linear_system_solver

    # (A + U*V')^-1
    function woodbury_identity(factored_A::abstract_linear_system_solver, U::AbstractArray, V::AbstractArray)
        invA_U = false
        capacitance_matrix = false
        try
          invA_U = ls_solve(factored_A, U)
          n = size(U,1);
          k = size(U,2);

          capacitance_matrix = (eye(k) + V * invA_U);
          factored_capacitance_matrix = lufact(capacitance_matrix);

          return new(U, V, invA_U, factored_capacitance_matrix, factored_A);
        catch e
            println("ERROR woodbury_identity")
            #@show full(factored_A._SparseMatrix)
            #@show capacitance_matrix
            #@show invA_U'
            throw(e)
        end
    end
end


function ls_solve(W::woodbury_identity, rhs::Array{Float64,1})
    my_temp = ls_solve(W.factored_A, rhs)
    return (my_temp - W.invA_U * (W.factored_capacitance_matrix \ (W.V * my_temp)))
end

function numerical_error(A, sol, rhs)
    return norm(A * sol - rhs,1);
end

function e_(i::Int64,n::Int64)
    @assert(i <= n)
    vec = zeros(n);
    vec[i] = 1.0;
    return vec;
end
