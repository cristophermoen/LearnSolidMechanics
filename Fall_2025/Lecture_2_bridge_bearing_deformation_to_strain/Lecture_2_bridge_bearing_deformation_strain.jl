using ForwardDiff

#Define deformation fields for a bridge bearing and calculate deformation gradients point by point.

H = 12.0 #in.
W = 27.0 #in.
μ = 0.5 #Poisson's ratio for natural rubber 

a = -2.0 #displacement at top of bearing, plane A-A-A-A
b = -μ * a  #maximum bulging or lateral displacement at midheight of bearing 

#point = point defined with coordinates [x, y, z]
                                        [1, 2, 3]

#define the vertical deformation field
v(point) = a * point[2] / H 

#define lateral deformation field 
u(point) = b * sin(π * point[2]/ H) * (1/(W/2)) * (point[1] - W/2)
w(point) = b * sin(π * point[2]/ H) * (1/(W/2)) * (point[3] - W/2)


#now calculate the gradient of the deformation fields at a point

#define the coordinates of the point
point = [W/2, 0.0, W]

#calculate deformation at this point, just to check things 
u(point)
v(point)
w(point)

#this gives ∂u/∂x, ∂u/∂y, ∂u/∂z at a point in the domain 
∇u_at_a_point = ForwardDiff.gradient(u, point)

#this gives ∂v/∂x, ∂v/∂y, ∂v/∂z at a point in the domain 
∇v_at_a_point = ForwardDiff.gradient(v, point)

#this gives ∂w/∂x, ∂w/∂y, ∂w/∂z at a point in the domain 
∇w_at_a_point = ForwardDiff.gradient(w, point)

#displacement gradient tensor
∇ = [∇u_at_a_point'
    ∇v_at_a_point'
    ∇w_at_a_point']


#strain tensor
ϵkl_at_a_point = 1/2 .* [∇[i,j] + ∇[j,i] for i=1:3, j=1:3]