#= License
Copyright 2019, 2020 (c) Yossi Bokor Katharine Turner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

__precompile__()


module DiscretePersistentHomologyTransform


#### Requirements ####
using CSV
using Hungarian
using DataFrames
using LinearAlgebra
using SparseArrays
using Eirene
using Plots

#### Exports ####

export 	PHT,
		Recenter,
		Direction_Filtration,
		Evaluate_Barcode,
		Total_Rank_Exact,
		Total_Rank_Grid,
		Total_Rank_Auto,
		Average_Rank_Grid,
		Create_Heat_Map,
		Set_Mean_Zero,
		Weighted_Inner_Product,
		Weighted_Inner_Product_Matrix,
		Principal_Component_Scores,
		Average_Discretised_Rank,
		Plot_Diagrams,
		Plot_Eirene_Diagram,
		unittest
		
		
#### First some functions to recenter the curves ####
function Find_Center(points)
	n_p = size(points,1)
	
	c_x = Float64(0)
	c_y = Float64(0)
	
	for i in 1:n_p
		c_x += points[i,1]
		c_y += points[i,2]
	end
	
	return Float64(c_x/n_p), Float64(c_y/n_p)
end #Find_Center

function Recenter(points; mode="centroid", directions=1, C=1)
	
	number_of_points = size(points,1)
	if mode == "centroid"
		points = convert(Array{Float64}, points)
		center = Find_Center(points)
	
		for i in 1:size(points)[1]
			points[i,1] = points[i,1] - center[1]
			points[i,2] = points[i,2] - center[2]
		end
		
		return points
	elseif mode == "center"
		dirs = Array{Float64}(undef, directions,2)
		for n in 1:directions
			dirs[n,1] = cos(n*pi/(directions/2))
			dirs[n,2] = sin(n*pi/(directions/2))
		end
		
		lambda = []
		
		for i in 1:directions
			heights = Array{Float64}(undef, 1, number_of_points)
			direction = dirs[i,:]
			for i in 1:number_of_points
				heights[i]= points[i,1]*direction[1] + points[i,2]*direction[2] #calculate heights in specificed direction
			end
			append!(lambda, minimum(heights))
		end
		
		K = directions/2 #need to check with Kate about this.
		
		ld =lambda.*dirs
		u = (1/K)*[sum(ld[:,1]), sum(ld[:,2])]
		
		points = points.-u'
		
		return points
	
	elseif mode =="scale"

		dirs = Array{Float64}(undef, directions,2)
		for n in 1:directions
			dirs[n,1] = cos(n*pi/(directions/2))
			dirs[n,2] = sin(n*pi/(directions/2))
		end
		
		lambda = []
		
		for i in 1:directions
			heights = Array{Float64}(undef, 1, number_of_points)
			direction = dirs[i,:]
			for i in 1:number_of_points
				heights[i]= points[i,1]*direction[1] + points[i,2]*direction[2] #calculate heights in specificed direction
			end
			append!(lambda, minimum(heights))
		end
		
		K = directions/2 #need to check with Kate about this.
		
		ld =lambda.*dirs
		u = (1/K)*[sum(ld[:,1]), sum(ld[:,2])]
		
		points = points.-u'
		
		L = -sum(lambda)
		
		sf = C/L
		
		return sf.*points
	
	
	end
	
end #Recenter


function Evaluate_Rank(barcode, point)
	
	n = size(barcode,1)
	count = 0
	
	if point[2] < point[1]
		#return 0
	else
		for i in 1:n
			if barcode[i,1] <= point[1]
				if barcode[i,2] >= point[2]
					count +=1
				end
			end
		end
		return count
	end
end

function Total_Rank_Exact(barcode)
	@assert size(barcode,2) == 2

	rks = []

	
	n = size(barcode,1)
	
	b = copy(barcode)
	reshape(b, 2,n)
	for i in 1:n
		for j in 1:n 
			if barcode[i,1] < barcode[j,1]
				if barcode[i,2] < barcode[j,2]
					if barcode[i,1] < barcode[j,2]
						b = vcat(b, [barcode[i,1] barcode[j,2]])
					end
				end
			end
		end
	end
	for i in 1:size(b,1)
		append!(rks, Evaluate_Rank(barcode, b[i,:]))
	end
	return b, rks

end


function Total_Rank_Grid(barcode, x_g, y_g) #the grid should be an array, with 0s in all entries below the second diagonal.
		
		@assert size(x_g) == size(y_g) #I should maybe change this as you don't REALLY need to use the same size grid.....

		n_g = size(x_g,1)
		rks = zeros(n_g,n_g)
		n_p = size(barcode,1)

		for i in 1:n_p
			point = barcode[i,:]
			x_i = findfirst(>=(point[1]), x_g)
			y_i = findfirst(<=(point[2]), y_g)
			for j in x_i:n_g-y_i+1
				for k in j:n_g-y_i+1
					rks[n_g-k+1,j] += 1
				end
			end
		end

	return rks

end

function Average_Rank_Grid(list_of_barcodes, x_g, y_g)

	
	rks = zeros(length(x_g),length(y_g))
	
	
	n_b = length(list_of_barcodes)
	
	for i in 1:n_b
		rk_i = Total_Rank_Grid(list_of_barcodes[i], x_g,y_g)
		rks = rks .+ rk_i
	end
	
	rks = rks/n_b
	
	return rks
end

function Average_Rank_Point(list_of_barcodes, x,y)
	
	rk = 0
	n_b = length(list_of_barcodes)
	if y >= x
	for i in 1:n_b
			rk += Evaluate_Rank(list_of_barcodes[i], [x,y])
		end
		return rk/n_b
	else
		 return 0
	end
end

function Create_Heat_Map(barcode, x_g, y_g)
	
	f(x,y) =  begin
					if x > y
						return 0
					else
						return Evaluate_Rank(barcode,[x,y])
					end
				end
				
	#Z = map(f, X, Y)

	p1 = contour(x_g, y_g, f, fill=true)

	return p1
end

#### Doing 0th Persistent Homology Using Union Find with using birth of the essential 1-cycle as the death for the first 0-cycle ####

function UF_PH(ord_pts::AbstractArray, direction::AbstractArray) #want ord_pts to be mx2 array
	n_p = size(ord_pts, 1)
	heights = Dict()
	bd = Dict{Int, Array}()
	parents = Dict{Int, Array}()
	for i in 1:n_p
		heights[i]=dot(ord_pts[i,:], direction)
	end
	appearance = sort(collect(heights), by=x->x[2])
	bd[appearance[1][1]]=[appearance[1][2]]
	parents[appearance[1][1]]= [appearance[1][1]]
	Seen = Set(appearance[1][1])
	for j in 2:n_p-1
		i = appearance[j][1]
		if i == n_p
			n_b = Set([i-1, 1])
		elseif i == 1
			n_b = Set([i+1, n_p])
		else
			n_b = Set([i-1, i+1])
		end
		seen = collect(intersect(n_b, Seen))
		if size(seen,1) == 0
			parents[i]=[i]
			bd[i] = [heights[i]]
		elseif size(seen,1) == 1
			for k in keys(parents)
				if size(collect(intersect(Set(parents[k]), Set(seen[1]))),1) !=0
					parents[k]=append!(parents[k], i)
				end
			end
		else
			pars = []
			for k in keys(parents)
				if size(collect(intersect(Set(parents[k]), Set(seen)))) != 0
					append!(pars, k)
				end
			end
			if heights[pars[1]] <= heights[pars[2]]
				parents[pars[1]] = append!(parents[pars[1]], parents[pars[2]])
				parents[pars[1]] = append!(parents[pars[1]], i)
				bd[pars[2]] = append!(bd[pars[2]], heights[i])
				delete!(parents, pars[2])
			else
				parents[pars[2]] = append!(parents[pars[2]], parents[pars[1]])
				parents[pars[2]] = append!(parents[pars[2]], [i])
				bd[pars[1]] = append!(bd[pars[1]], heights[i])
				delete!(parents, pars[1])
			end
		end
		push!(Seen, i)
	end
	@assert size(collect(Seen),1) == n_p-1
	i = appearance[n_p][1]
	@assert size(collect(keys(parents)),1) == 1
	bd[collect(keys(parents))[1]] = append!(bd[collect(keys(parents))[1]], appearance[n_p][2])
	push!(Seen, i)
	@assert size(collect(Seen),1) == n_p
	pers_pts = collect(keys(bd))
	for i in pers_pts
		if bd[i][1] == bd[i][2]
			delete!(bd, i)
		end
	end
	return bd
end


# Let us do PCA for the rank functions using Kate and Vanessa's paper.
# So, I first need to calculate the pointwise norm

function Set_Mean_Zero(discretised_ranks)
	n_r = length(discretised_ranks)
	grid_size = size(discretised_ranks[1])
	
	for i in 1:n_r
		@assert size(discretised_ranks[i]) == grid_size
	end
	
	mu = zeros(grid_size)
	
	for i in 1:n_r
		mu = mu .+ discretised_ranks[i]
	end
	
	mu = mu./n_r
	
	normalised = []
	for i in 1:n_r
		append!(normalised, [discretised_ranks[i] .- mu])
	end
	return normalised
end

function Weighted_Inner_Product(disc_rank_1, disc_rank_2, weights)

	wip = sum((disc_rank_1.*disc_rank_2).*weights)

	return wip
end

function Weighted_Inner_Product_Matrix(discretised_ranks, weights)
	n_r = length(discretised_ranks)
	D = Array{Float64}(undef, n_r, n_r)
	
	for i in 1:n_r
		for j in i:n_r
			wip = Weighted_Inner_Product(discretised_ranks[i], discretised_ranks[j], weights)
			D[i,j] = wip
			D[j,i] = wip
		end
	end
	
	return D
end

function Principal_Component_Scores(inner_prod_matrix, dimension)
	F = LinearAlgebra.eigen(inner_prod_matrix, permute = false, scale=false) # this sorts the eigenvectors in ascending order
	n_r = size(inner_prod_matrix,1)
	lambda = Array{Float64}(undef, 1,dimension)
	w = Array{Float64}(undef, size(F.vectors)[1],dimension)
	n_v = length(F.values)

	for i in 1:dimension
		lambda[i] = F.values[n_v-i+1]
		w[:,i] = F.vectors[:,n_v-i+1]
	end
	
	s = Array{Float64}(undef, n_r,dimension)
	
	for i in 1:size(inner_prod_matrix,1)
		for j in 1:dimension
			
			den = sqrt(sum([w[k,j]*sum(w[l,j]*inner_prod_matrix[k,l] for l in 1:n_r) for k in 1:n_r]))
			numerator = sum(w[:,j].*inner_prod_matrix[:,i])
			s[i,j] = numerator/den
		end
	end
	return s
end

function Average_Discretised_Rank(list_of_disc_ranks)
	average = Array{Float64}(undef, size(list_of_disc_ranks[1]))
	n_r = length(list_of_disc_ranks)
	
	for i in n_r
		average = average .+ list_of_disc_ranks[i]
	end
	
	return average/n_r
end


function Direction_Filtration(ordered_points, direction; out = "barcode", one_cycle = false )
	number_of_points = length(ordered_points[:,1]) #number of points
	heights = zeros(number_of_points) #empty array to be changed to heights for filtration
	fv = zeros(2*number_of_points) #blank fv Eirene
	for i in 1:number_of_points
		heights[i]= ordered_points[i,1]*direction[1] + ordered_points[i,2]*direction[2] #calculate heights in specificed direction
	end
	
	for i in 1:number_of_points
		fv[i]= heights[i] # for a point the filtration step is the height
	end
	
	for i in 1:(number_of_points-1)
		fv[(i+number_of_points)]=maximum([heights[i], heights[i+1]]) # for an edge between two adjacent points it enters when the 2nd of the two points does
	end
	
	fv[2*number_of_points] = maximum([heights[1] , heights[number_of_points]]) #last one is a special snowflake
	dv = [] # template dv for Eirene
	
	for i in 1:number_of_points
		append!(dv,0) # every point is 0 dimensional
	end
	
	for i in (1+number_of_points):(2*number_of_points)
		append!(dv,1) # edges are 1 dimensional
	end
	
	D = zeros((2*number_of_points, 2*number_of_points))
	
	for i in 1:number_of_points
		D[i,(i+number_of_points)]=1 # create boundary matrix and put in entries
	end
	
	for i in 2:(number_of_points)
		D[i, (i+number_of_points-1)]=1 # put in entries for boundary matrix
	end
	
	D[1, (2*number_of_points)]=1
	
	ev = [number_of_points, number_of_points] # template ev for Eirene
	
	S  = sparse(D) 	# converting as required for Eirene
	rv = S.rowval 	# converting as required for Eirene
	cp = S.colptr	# converting as required for Eirene
	C = Eirene.eirene(rv=rv,cp=cp,ev=ev,fv=fv) # put it all into Eirene
	
	if out == "barcode"
		if one_cycle == true
			return barcode(C, dim=0), maximum(heights)
		else
			return barcode(C, dim=0)
		end
	else
		if one_cycle == true
			return C, maximum(heights)
		else
			return C
		end
	end
end #Direction_Filtration

# compare multiple persistence diagrams on the same plot)
function Plot_Diagrams(diagrams::Array)
	n_d = length(diagrams)
	pyplot()
	
	scatter(diagrams[1][:,1], diagrams[1][:,2], label="Diagram 1")
	if n_d != 1
		for i in 2:n_d-1
			scatter!(diagrams[i][:,1], diagrams[i][:,2], label="Diagram $i")
		end
	end
	scatter!(diagrams[n_d][:,1], diagrams[n_d][:,2], label="Diagram $n_d")
end

function Plot_Diagrams(diagrams::Array, labels::Array{String})
	n_d = length(diagrams)
	pyplot()
	
	@assert n_d == length(labels)
	scatter(diagrams[1][:,1], diagrams[1][:,2], label=labels[1])
	if n_d != 1
		for i in 2:n_d-1
			scatter!(diagrams[i][:,1], diagrams[i][:,2], label=labels[i])
		end
	end
	scatter!(diagrams[n_d][:,1], diagrams[n_d][:,2], label=labels[n_d])
end


function Plot_Eirene_Diagram(C)
	Eirene.plotpersistencediagram_pjs(C, dim=0)
end


#### Wrapper for the PHT function ####

function PHT(curve_points, directions;method="Eirene", out="barcode", one_cycle = false) #accepts an ARRAY of points
	
	if typeof(directions) ==  Int64
		println("auto generating directions")
		dirs = Array{Float64}(undef, directions,2)
		for n in 1:directions
			dirs[n,1] = cos(n*pi/(directions/2))
			dirs[n,2] = sin(n*pi/(directions/2))
		end
		println("Directions are:", dirs)
	else
		println("using directions provided")
		dirs = copy(directions)
	end
	pht = []
	
	if one_cycle == true
		c_1 = []
	end
	
	for i in 1:size(dirs,1)
		if method == "Eirene"
			if one_cycle == true
				pd,c_1 = Direction_Filtration(curve_points, dirs[i,:], one_Cycle = true)
				append!(cycle_1, c_1)
				pht = vcat(pht, [pd])
			else
				pd = Direction_Filtration(curve_points, dirs[i,:])
				pht = vcat(pht, [pd])
			end
		elseif method == "UnionFind"
			bd = UF_PH(curve_points, dirs[i,:])
			pd=Array{Float64}(undef, 0,2)
			for k in values(bd)
				pd = vcat(pd, k')
			end
			pht=vcat(pht, [pd])
		else
			println("Error: Pleasure ensure you specify either Eirene or UnionFind as your method.")
			return
		end
	end
	if one_cycle == true
		return pht, cycle_1
		
	else
		return pht
	end
end


#### Wrapper for PCA ####
function PCA(ranks, dimension, weights)
	
	normalised = Set_Mean_Zero(ranks)

	D = Weighted_InnerProd_Matrix(normalised, weights)
	
	return Principal_Component_Scores(D, dimension)
end

#### Unittests ####
function test_1()
	return PHT([0,0,0],0)
end

function test_2()
	pht = PHT([1 1; 5 5], 1)
	if pht == [0.9999999999999998 Inf]
		return []
	else
		println("Error: test_2, pht = ")
		return pht
	end
end

function unittest()

	x = Array{Any}(undef, 2)
	
	x[1] = test_1()
	x[2] = test_2()
	
	for p 	= 	1:length(x)
		if !isempty(x[p])
			println(p)
			return x
		end
	end
	return []
end

end# module
