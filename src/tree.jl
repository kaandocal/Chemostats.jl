using ArgCheck
using UnPack

const OffspringType{T} = Union{Nothing, Tuple{T}, Tuple{T, T}}

struct AncestryTree{T}
    parents::WeakKeyDict{T,T}
    children::WeakKeyDict{T,OffspringType{T}}
    leaves::Vector{T}
    save_ancestors::Bool 
    save_children::Bool
    save_leaves::Bool
end

function PopTree{T}(; save_ancestors=false, save_children=false, save_leaves=false) where T
    PopTree(WeakKeyDict{T,T}(), WeakKeyDict{T,OffspringType{T}}(), T[], save_ancestors, save_children, save_leaves)
end 

function parent(tree::PopTree{T}, obj::T) where {T}
    get(tree.parents, obj, missing)
end 

function set_parent!(tree::PopTree{T}, obj::T, parent::T) where {T}
    @check !haskey(tree.parents, obj) "Attempting to assign parent to cell which already has parent"

    tree.parents[obj] = parent
end 

function children(tree::PopTree{T}, obj::T) where {T} 
    get(tree.children, obj, missing)
end 

function add_offspring!(tree::PopTree{T}, parent::T, children::OffspringType{T}) where {T}
    @check !haskey(tree.children, parent) "Attempting to assign children to cell which already has children"

    if tree.save_ancestors
        for cell in children
            set_parent!(tree, cell, parent)
        end 
    end 

    if tree.save_children
        tree.children[parent] = children
    end

    nothing
end 

function add_leaf!(tree::PopTree{T}, obj::T) where {T}
    @check !haskey(tree.children, obj) "Attempting to assign children to cell which already has children"

    tree.children[obj] = nothing 

    if tree.save_leaves 
        push!(tree.leaves, obj)
    end 
end

# function clone!(tree, obj)
# end 

# Need lifetime support
# function snapshot(tree, t)
# end 

# function alive_at(node, t)


struct BackwardsIterator{T}
    tree::PopTree{T}
    obj::T
end 

ancestors(tree::PopTree{T}, obj::T) where T = BackwardsIterator(tree, obj)

function Base.iterate(iter::BackwardsIterator{T}, obj::T=iter.obj) where T
    @unpack tree = iter 

    if haskey(tree.parents, obj)
        parent = tree.parents[obj]
        parent, parent
    else 
        nothing 
    end 
end 
