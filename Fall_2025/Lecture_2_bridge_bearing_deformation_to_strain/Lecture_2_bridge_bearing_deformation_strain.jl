using ForwardDiff

#Define deformation fields for a bridge bearing and calculate deformation gradients point by point.

H = 12.0 #in.
W = 27.0 #in.
μ = 0.5 #Poisson's ratio for natural rubber 

a = -1.0 #displacement at top of bearing, plane A-A-A-A
b = -μ * a  #maximum bulging or lateral displacement at midheight of bearing 

#point = point defined with coordinatees [x, y, z]


#define the vertical deformation field
v(point) = a * point[2] / H 

#define lateral deformation field 
u(point) = b * sin(π * point[2]/ H) * (1/(W/2)) * (point[1] - W/2)


#now calculate the gradient of the deformation fields at a point

#define the coordinates of the point
point = [W, H/2, 0.0]

#calculate u at this point, just to check things 
u(point)

#calculate v at this point, just to check things 
v(point)

#this gives ∂u/∂x, ∂u/∂y, ∂u/∂z at a point in the domain 
∇u_at_a_point = ForwardDiff.gradient(u, point)

#this gives ∂v/∂x, ∂v/∂y, ∂v/∂z at a point in the domain 
∇v_at_a_point = ForwardDiff.gradient(v, point)



