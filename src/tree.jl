using ArgCheck
using UnPack

const OffspringType{T} = Union{Nothing, Tuple{T}, Tuple{T, T}}

mutable struct PopTree{T}
    parents::Dict{T,T}
    leaves::Vector{T}
    save_ancestors::Bool
    save_leaves::Bool
    last_sweep_size::Int
end

function PopTree{T}(; save_ancestors=false, save_leaves=false) where T
    PopTree{T}(Dict{T,T}(), T[], save_ancestors, save_leaves, 0)
end

function parent(tree::PopTree{T}, obj::T) where {T}
    get(tree.parents, obj, missing)
end 

function set_parent!(tree::PopTree{T}, obj::T, parent::T) where {T}
    @check !haskey(tree.parents, obj) "Attempting to assign parent to cell which already has parent"

    tree.parents[obj] = parent
end 

function add_offspring!(tree::PopTree{T}, parent::T, children::Union{Nothing,Tuple{Vararg{T}}}) where {T}
    if tree.save_ancestors
        for cell in children
            set_parent!(tree, cell, parent)
        end
    end

    nothing
end 

function add_leaf!(tree::PopTree{T}, obj::T) where {T}
    if tree.save_leaves
        push!(tree.leaves, obj)
    end
end

function simplify!(tree::PopTree, samples)
    if tree.save_ancestors && length(tree.parents) > 2 * tree.last_sweep_size
        prune!(tree, samples)
    end
end

"""
    prune!(tree::PopTree, samples)

Prunes `tree.parents` to eliminate dead leaves.
"""
function prune!(tree::PopTree{T}, samples) where T
    new_parents = Dict{T,T}()
    sizehint!(new_parents, tree.last_sweep_size)

    for s in samples
        c = s
        while haskey(tree.parents, c) && !haskey(new_parents, c)
            p = tree.parents[c]
            new_parents[c] = p
            c = p
        end
    end

    tree.parents = new_parents
    tree.last_sweep_size = length(new_parents)
    tree
end

struct BackwardsIterator{T}
    tree::PopTree{T}
    obj::T
end 

ancestors(tree::PopTree{T}, obj::T) where T = BackwardsIterator(tree, obj)
Base.eltype(::BackwardsIterator{T}) where T = T
Base.IteratorSize(::BackwardsIterator) = Base.SizeUnknown()

function Base.iterate(iter::BackwardsIterator{T}, obj::T=iter.obj) where T
    @unpack tree = iter 

    if haskey(tree.parents, obj)
        parent = tree.parents[obj]
        parent, parent
    else 
        nothing 
    end 
end 
