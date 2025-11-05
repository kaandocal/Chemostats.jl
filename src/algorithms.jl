# Adapted from DataStructures.jl 
_getindex(queue::BinaryHeap, i::Integer) = queue.valtree[i]

function _force_up!(queue::BinaryHeap, i::Integer)
    x = queue.valtree[i]

    @inbounds while i > 1
        j = DataStructures.heapparent(i)
        queue.valtree[i] = queue.valtree[j]
        i = j
    end

    queue.valtree[i] = x
end 

function _popat!(queue::BinaryHeap, idx::Integer)
    _force_up!(queue, idx)
    pop!(queue)
end 


abstract type AbstractAlgorithm end 

struct Forward <: AbstractAlgorithm
end 

get_δ(::Forward) = 0.
init!(alg::Forward, int) = nothing
update_algorithm!(::Forward, int) = nothing

function add_offspring!(queue, cells, ::Forward)
    isempty(cells) && return 0
    push!(queue, first(cells))
    return 1
end

struct Thin <: AbstractAlgorithm
    δ::Float64 
end 

get_δ(alg::Thin) = alg.δ
init!(alg::Thin, int) = nothing
update_algorithm!(::Thin, int) = nothing

function add_offspring!(queue, cells, ::Thin)
    for cell in cells
        push!(queue, cell)
    end

    return length(cells)
end

struct Strict <: AbstractAlgorithm
    L::Int
end 

Strict(cells::AbstractVector) = Strict(length(cells))
get_δ(::Strict) = 0.
init!(alg::Strict, int) = update_algorithm!(alg, int)

function update_algorithm!(alg::Strict, int) 
    if isempty(int.queue)
        @warn "No cells left in strict chemostat"
        int.retcode = ReturnCode.Unstable
        return
    end 

    # Ensure we have exactly L clones 
    N_start = length(int.queue)

    while length(int.queue) < alg.L
        duplicate_random!(int.queue, int.t)
    end

    while length(int.queue) > alg.L
        j = sample(1:length(int.queue))
        _popat!(int.queue, j)
    end 

    int.log_f += log(N_start) - log(alg.L)
end 

function duplicate_random!(queue::BinaryHeap, t)
    i = sample(1:length(queue))
    push!(queue, duplicate_cell(queue.valtree[i], t))
end 

function add_offspring!(queue, cells, alg::Strict)
    @check length(queue) == alg.L - 1

    nsim = 0

    for cell in cells
        push!(queue, cell)
        nsim += 1
    end 

    while length(queue) > alg.L
        j = sample(1:length(queue))
        if _getindex(queue, j) in cells
            nsim -= 1
        end 

        _popat!(queue, j)
    end 

    return nsim
end

# struct Lax?