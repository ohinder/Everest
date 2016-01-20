# QCQP


#using Everest

include("src/loadcode.jl");
test_problem_nlo("BRANDY")

##############################################################
#############################################################
testmat = rand(2,2);
testmat[1,1] = 1e8;
testmat[1,2] = 1e8;
testmat[2,1] = 1e8 + 1;
testmat[2,2] = 1e8
testrhs = rand(2)
sol1 = testmat \ testrhs
norm(testmat*sol1-testrhs);

sol2 = lufact(testmat) \ testrhs
norm(testmat*sol2-testrhs);
norm(testmat*inv(testmat)*testrhs-testrhs);


function test(R::SparseMatrixCSC{Float64,Int64},v::Array{Float64,1})
  println("cholfac")
  @time F = cholfact(R + 1e-7 * speye(size(R,1)));

  println("solve")

  @time begin
    for i = 1:100
        sol = F \ v
    end
  end
end

function test2(A::SparseMatrixCSC{Float64,Int64},v::Array{Float64,1})
  println("lufac")
  @time F = lufact([speye(size(A,2))]);

  println("solve")

  @time begin
    for i = 1:50
        F \ v # compare with mumps
    end
  end
end

function test3(R::SparseMatrixCSC{Float64,Int64},V::AbstractMatrix)
  println("cholfac")
  @time F = cholfact(R + 1e-7 * speye(size(R,1)));

  println("solve")


  @time sol = F \ V
end

A,b,c = get_netlib_problem("QAP15")
@time R = A*A';

v = ones(size(R,1));
test(R,v)

siz = 100;
rhs = zeros(size(R,1),siz);
rhs[1,:] = ones(siz)
test3(R,rhs)

using MUMPS
