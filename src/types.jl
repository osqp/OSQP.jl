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

# Internal C type for info
# N.B. This is not the one returned by the user!
struct CInfo
	iter::Clong
	# We need to allocate 32 bytes for a character string, so we allocate 256 bits
	# of integer instead
	# TODO: Find a better way to do this
	# status1::Int128
	# status2::Int128
	status::NTuple{32, Cchar}
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


end

function Settings()
	s = Ref{OSQP.Settings}()
	ccall((:set_default_settings, OSQP.osqp), Void,
	      (Ref{OSQP.Settings},), s)
	return s[]
end

# function Settings(settings_dict::Dict{Symbol, Any})
function Settings(settings::Array{Any, 1})
	default_settings = OSQP.Settings()

	settings_dict = Dict{Symbol, Any}()
	if !isempty(settings)
		for (key, value) in settings
			settings_dict[key] = value
		end
	end

	# Get dictionary with elements of detauls and user settings
	settings_list = [setting in keys(settings_dict) ?
			 convert(fieldtype(typeof(default_settings), setting), settings_dict[setting]) :
			 getfield(default_settings, setting)
			 for setting in fieldnames(default_settings)]

	# Create new settings with new dictionary
	s = OSQP.Settings(settings_list...)
	return s
	
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
	info::Ptr{OSQP.CInfo}

	timer::Ptr{Void}
	first_run::Clong
	summary_printed::Clong

end


struct Info
	iter::Int64
	status::Symbol
	status_val::Int
	status_polish::Int
	obj_val::Float64
	pri_res::Float64
	dua_res::Float64
	setup_time::Float64
	solve_time::Float64
	polish_time::Float64
	run_time::Float64

	function Info(cinfo::CInfo) 
		status = OSQP.status_map[cinfo.status_val]
		return new(cinfo.iter, status,
			   cinfo.status_val,
			   cinfo.status_polish,
			   cinfo.obj_val,
			   cinfo.pri_res,
			   cinfo.dua_res,
			   cinfo.setup_time,
			   cinfo.solve_time,
			   cinfo.polish_time,
			   cinfo.run_time)
	end
end

struct Results
	x::Vector{Float64}
	y::Vector{Float64}
	info::OSQP.Info
end
