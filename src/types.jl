# Types defined in types.h
# https://github.com/oxfordcontrol/osqp/blob/master/include/types.h



# export Ccsc, OSQPData, OSQPSettings  # TODO: remove


struct Ccsc
	nzmax::Clong
	m::Clong
	n::Clong
	p::Ptr{Clong}
	i::Ptr{Clong}
	x::Ptr{Cdouble}
	nz::Clong


	# Constructor from SparseMatrixCSC
	function Ccsc(M::SparseMatrixCSC)

		# Get dimensions
		m = M.m
		n = M.n

		# Get vectors of data, rows indices and column pointers
		x = convert(Array{Float64, 1}, M.nzval)
		# C is 0 indexed
		i = convert(Array{Clong, 1}, M.rowval - 1)  
		# C is 0 indexed
		p = convert(Array{Clong, 1}, M.colptr - 1)  

		new(length(M.nzval), m, n, pointer(p), pointer(i), pointer(x), -1)
	end
end

struct Solution
	x::Ptr{Cdouble}
	y::Ptr{Cdouble}
end

struct Info
	iter::Clong
	# We need to allocate 32 bytes for a character string, so we allocate 256 bits
	# of integer instead
	# TODO: Find a better way to do this
	status1::Int128
	status2::Int128
	status_val::Clong
	status_polish::Clong
	obj_val::Cdouble
	pri_res::Cdouble
	dua_res::Cdouble
	setup_time::Cdouble
	solve_time::Cdouble
	polish_time::Cdouble
	run_time::Cdouble
end

struct Data
	n::Clong
	m::Clong
	P::Ptr{Ccsc}
	A::Ptr{Ccsc}
	q::Ptr{Cdouble}
	l::Ptr{Cdouble}
	u::Ptr{Cdouble}

	function Data(n::Int,
			  m::Int,
			  P::SparseMatrixCSC,
			  q::Vector{Float64}, 
			  A::SparseMatrixCSC, 
			  l::Vector{Float64},
			  u::Vector{Float64})
	
		Pcsc = OSQP.Ccsc(P)
		Acsc = OSQP.Ccsc(A)
		new(n, m, pointer([Pcsc]), pointer([Acsc]), 
		    pointer(q), pointer(l), pointer(u))
	end
end


struct Settings
	rho::Cdouble
	sigma::Cdouble
	scaling::Clong
	scaling_iter::Clong
	max_iter::Clong
	eps_abs::Cdouble
	eps_rel::Cdouble
	eps_prim_inf::Cdouble
	eps_dual_inf::Cdouble
	alpha::Cdouble
	delta::Cdouble
	polish::Clong
	pol_refine_iter::Clong
	verbose::Clong
	auto_rho::Clong
	early_terminate::Clong
	early_terminate_interval::Clong
	warm_start::Clong


	function Settings()
		s = Ref{OSQP.Settings}() 
	    ccall((:set_default_settings, OSQP.osqp), Void, 
		  (Ref{OSQP.Settings},), s)

	    # Ref Settings can be accessed using stgs = s[]
	    return s[]
	end
end


struct Workspace
	data::Ptr{OSQP.Data}
	priv::Ptr{Void}
	pol::Ptr{Void}

	x::Ptr{Cdouble}
	y::Ptr{Cdouble}
	z::Ptr{Cdouble}
	xz_tilde::Ptr{Cdouble}
	x_prev::Ptr{Cdouble}
	z_prev::Ptr{Cdouble}
	delta_y::Ptr{Cdouble}
	Atdelta_y::Ptr{Cdouble}
	delta_x::Ptr{Cdouble}
	Pdelta_x::Ptr{Cdouble}
	Adelta_x::Ptr{Cdouble}
	D_temp::Ptr{Cdouble}
	D_temp_A::Ptr{Cdouble}
	E_temp::Ptr{Cdouble}

	settings::Ptr{OSQP.Settings}
	scaling::Ptr{Void}
	solution::Ptr{OSQP.Solution}
	info::Ptr{OSQP.Info}

	timer::Ptr{Void}
	first_run::Clong
	summary_printed::Clong

end


struct Results
	x::Vector{Float64}
	y::Vector{Float64}
	# info::OSQP.Info
end
