using ArgCheck
using UnPack

# Cell -> user-defined data
# Each cell has a unique ID 

const OffspringType{T} = Union{Nothing, Tuple{T}, Tuple{T, T}}

struct PopTree{T,ID}
    data::Dict{ID,T}
    parents::Dict{ID,ID}
    children::Dict{ID,OffspringType{ID}}
    leaves::Vector{ID}
    save_lineages::Bool 
    save_leaves::Bool
end

function PopTree(roots::AbstractVector; save_lineages=false, save_leaves=false)
    data = Dict(objectid(root) => root for root in roots)
    ID = eltype(keys(data))
    PopTree(data, Dict{ID,ID}(), Dict{ID,OffspringType{ID}}(), ID[], save_lineages, save_leaves)
end 

function add_cell!(tree, cell) 
    id = objectid(cell)
    @check !haskey(tree.data, id) "Adding duplicate cell to population tree"
    tree.data[id] = cell
    tree
end 

function getid(tree::PopTree{T}, obj::T) where T 
    id = objectid(obj)
    haskey(tree.data, id) || raise(BoundsError("Attempting to index cell not in tree"))
    id 
end 

function parent(tree::PopTree{T}, obj::T) where {T}
    pid = get(tree.parents, getid(tree, obj), missing)

    ismissing(pid) && return missing 
    tree.data[pid]
end 

function set_parent!(tree::PopTree{T,ID}, obj::T, pid::ID) where {T,ID}
    id = getid(tree, obj)
    @check !haskey(tree.parents, id) "Attempting to assign parent to cell which already has parent"

    tree.parents[id] = pid
end 

function children(tree::PopTree{T}, obj::T) where {T} 
    ch = get(tree.children, getid(tree, obj), missing)

    (isnothing(ch) || ismissing(ch)) && return ch
    map(id -> tree.data[id], ch)
end 

function add_offspring!(tree::PopTree{T}, parent::T, children::OffspringType{T}) where {T}
    pid = getid(tree, parent)
    @check !haskey(tree.children, pid) "Attempting to assign children to cell which already has children"

    if isnothing(children)
        tree.children[pid] = nothing 
    else 
        for cell in children
            add_cell!(tree, cell)
            tree.save_lineages && set_parent!(tree, cell, pid)
        end 
        
        tree.children[pid] = map(cell -> getid(tree, cell), children)
    end

    nothing
end 


function die!(tree::PopTree{T}, obj::T) where {T}
    id = getid(tree, obj)
    @check !haskey(tree.children, id) "Attempting to assign children to cell which already has children"

    tree.children[id] = nothing 

    if tree.save_leaves 
        push!(tree.leaves, getid(tree, obj))
    end 
end

# function clone!(tree, obj)
# end 

# Need lifetime support
# function snapshot(tree, t)
# end 

# function alive_at(node, t)


struct BackwardsIterator{ID,TT <: PopTree{T,ID} where {T}}
    tree::TT
    id::ID 
end 

ancestors(tree::PopTree{T}, obj::T) where T = BackwardsIterator(tree, getid(tree, obj))

function Base.iterate(iter::BackwardsIterator{ID}, id::ID=iter.id) where {ID}
    @unpack tree = iter 

    if haskey(tree.parents, id)
        pid = tree.parents[id]
        tree.data[pid], pid
    else 
        nothing 
    end 
end 
