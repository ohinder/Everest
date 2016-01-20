# TO DO:
# - test final values
# - using fact check

# test the linear program solver on a simple problem

function trivial_lp(settings::class_settings)
	# trivial_test(A, b, c,  correct_status, problem_name, verbose)
	trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], :locally_optimal, "LP-OPTIMAL-A", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [1.0, 5.0], :locally_optimal, "LP-OPTIMAL-B", settings);
	trivial_test(sparse([1.0 5.0]), [10.0], [1.0, -1.0], :locally_optimal, "LP-OPTIMAL-C", settings);
	trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], :locally_optimal, "LP-OPTIMAL-D", settings);

	trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], :locally_primal_infeasible, "LP-INFEASIBLE-A", settings);
	trivial_test(sparse([1.0 5.0]), [-2.0], [-1.0, 1.0], :locally_primal_infeasible, "LP-INFEASIBLE-B", settings);
	trivial_test(sparse([1.0 1.0]), [-5.0], [0.0, 0.0], :locally_primal_infeasible, "LP-INFEASIBLE-C", settings);
	trivial_test(sparse([1.0 1.0]), [-0.001], [3.0, 1.0], :locally_primal_infeasible, "LP-INFEASIBLE-D", settings);

	trivial_test(spzeros(0,1), EMPTY_ARRAY, [-1.0], :locally_dual_infeasible, "LP-UNBOUNDED-A", settings);
	trivial_test(sparse([1.0 -5.0]), [0.0], [0.0, -1.0], :locally_dual_infeasible, "LP-UNBOUNDED-B", settings);
	trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], :locally_dual_infeasible, "LP-UNBOUNDED-C", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [0.0, -0.01], :locally_dual_infeasible, "LP-UNBOUNDED-D", settings);
end


function trivial_qp(settings::class_settings)
	Q1 = sparse([ [1.0 -1.0]; [-1.0 1.0] ]);
	Q2 = sparse(ones(2,2));
	Q3 = sparse([ [1.0 0.0]; [0.0 0.0] ]);

	# trivial_test(A, b, c, Q,  correct_status, problem_name, verbose)
	trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], speye(2), :locally_optimal, "QP-OPTIMAL-A", settings);
	trivial_test(sparse([1.0 -1.0]), [0.0], [-1.0, -1.0], Q3, :locally_optimal, "QP-OPTIMAL-B", settings);
	trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], Q2, :locally_optimal, "QP-OPTIMAL-C", settings);
	trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], Q1, :locally_optimal, "QP-OPTIMAL-D", settings);

	trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], speye(2), :locally_primal_infeasible, "QP-INFEASIBLE-A", settings);
	trivial_test(sparse([1.0 5.0]), [-1.0], [-1.0, -1.0], Q3, :locally_primal_infeasible, "QP-INFEASIBLE-B", settings);
	trivial_test(sparse([1.0 3.0]), [-1.0], [0.5, -2.0], Q2, :locally_primal_infeasible, "QP-INFEASIBLE-C", settings);
	trivial_test(sparse([1.0 1.0]), [-0.0001], [1.0, -1.0], Q1, :locally_primal_infeasible, "QP-INFEASIBLE-D", settings);


	trivial_test(sparse([0.0 0.0]), [0.0], [0.0, -1.0], Q3, :locally_dual_infeasible, "QP-UNBOUNDED-A", settings, true);
	trivial_test(sparse([1.0 0.0]), [5.0], [0.0, -1.0], Q3, :locally_dual_infeasible, "QP-UNBOUNDED-B", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [0.5, -2.0], Q1, :locally_dual_infeasible, "QP-UNBOUNDED-C", settings);
	trivial_test(sparse([1.0 -1.0]), [0.0], [0.0, -0.01], Q1, :locally_dual_infeasible, "QP-UNBOUNDED-D", settings);
end

function trivial_ncqp(settings::class_settings)
	Q1 = sparse([ [-1.0 1.0]; [1.0 -1.0] ]);
	Q2 = sparse(ones(2,2));
	Q3 = sparse([ [0.0 1.0]; [1.0 0.0] ]);


	# trivial_test(A, b, c, Q,  correct_status, problem_name, verbose)
	if false
		trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], speye(2), :locally_optimal, "NCQP-OPTIMAL-A", settings);
		trivial_test(sparse([1.0 1.0]), [1.0], [-1.0, -1.0], Q3, :locally_optimal, "NCQP-OPTIMAL-B", settings);
		trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], Q2, :locally_optimal, "NCQP-OPTIMAL-C", settings);
		trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], Q1, :locally_optimal, "NCQP-OPTIMAL-D", settings);

		trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], speye(2), :locally_primal_infeasible, "NCQP-INFEASIBLE-A", settings);
		trivial_test(sparse([1.0 5.0]), [-1.0], [-1.0, -1.0], Q3, :locally_primal_infeasible, "NCQP-INFEASIBLE-B", settings);
		trivial_test(sparse([1.0 3.0]), [-1.0], [0.5, -2.0], Q2, :locally_primal_infeasible, "NCQP-INFEASIBLE-C", settings);
		trivial_test(sparse([1.0 1.0]), [-0.0001], [1.0, -1.0], Q1, :locally_primal_infeasible, "NCQP-INFEASIBLE-D", settings);


		trivial_test(sparse([0.0 0.0]), [0.0], [0.0, -1.0], Q3, :locally_dual_infeasible, "NCQP-UNBOUNDED-A", settings);
		trivial_test(sparse([1.0 0.0]), [5.0], [0.0, -1.0], Q3, :locally_dual_infeasible, "NCQP-UNBOUNDED-B", settings);
		trivial_test(sparse([1.0 -1.0]), [1.0], [0.5, -2.0], Q1, :locally_dual_infeasible, "NCQP-UNBOUNDED-C", settings);
	  trivial_test(sparse([1.0 -1.0]), [1.0], [0.0, -0.01], Q1, :locally_dual_infeasible, "NCQP-UNBOUNDED-D", settings, true);
	end

end

 include("../src/loadcode.jl");
settings = class_settings();
settings.max_it = 100;
trivial_lp(settings);
trivial_qp(settings);




#trivial_lp(settings);
#trivial_qp(settings);
