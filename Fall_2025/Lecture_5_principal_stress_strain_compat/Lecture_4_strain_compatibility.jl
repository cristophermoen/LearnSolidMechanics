using ForwardDiff, LinearAlgebra

#Define deformation fields for a bridge bearing and calculate deformation gradients point by point.

H = 12.0 #in.
W = 27.0 #in.
ν = 0.4999 #Poisson's ratio for natural rubber 
E = 275.0 #psi 
G = E / (2 * (1 + ν)) #psi 

a = -0.5 #displacement at top of bearing, plane A-A-A-A
b = -ν * a  #maximum bulging or lateral displacement at midheight of bearing 

#point = point defined with coordinates [x, y, z]
                                        [1, 2, 3]

#define the vertical deformation field
v(point) = a * point[2] / H 

#define lateral deformation fields
u(point) = b * sin(π * point[2]/ H) * (1/(W/2)) * (point[1] - W/2)
w(point) = b * sin(π * point[2]/ H) * (1/(W/2)) * (point[3] - W/2)


#now calculate the gradient of the deformation fields at a point

#define the coordinates of the point
point = [W/2, 0.0, W]


###### EXAMPLE!!
# Define your function
f(x, y) = x^4 * y^3
x_val = 1.0
y_val = 1.0

# Compute the mixed derivative ∂²f/∂y∂x
mixed_derivative = ForwardDiff.derivative(y -> ForwardDiff.derivative(x -> f(x, y), x_val), y_val)
println(mixed_derivative) # Output: 12.0

mixed_derivative = ForwardDiff.derivative(y -> ForwardDiff.derivative(x -> f(x, y), x_val), y_val)

#####

u(x, y) = b * sin(π * y/ H) * (1/(W/2)) * (x - W/2)

#∂2ex/∂y2 = ∂2/∂y2(∂u/∂x)
mixed_derivative = ForwardDiff.derivative(y -> ForwardDiff.derivative(y -> ForwardDiff.derivative(x -> u(x, y), point[1]), point[2]), point[2])