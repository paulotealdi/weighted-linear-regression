from math import sqrt

x_list = [1.900,2.400,3.000,3.600,4.200]
y_list = [2.100,2.700,3.400,4.100,4.800]
ux_list = [0.058, 0.058, 0.058, 0.058, 0.058]
uy_list = [0.058, 0.058, 0.058, 0.058, 0.058]

def mult(arrA, arrB):
    arrMult = []
    for i in range(len(arrA)):
        arrMult.append(arrA[i]*arrB[i])
    return arrMult

def sum(arr):
    sum = 0
    for val in arr:
        sum += val
    return sum

def square(list):
    return [number ** 2 for number in list]

def inv_square(list):
    return [(1/(number ** 2)) for number in list]

def reglin_conv(x_list, y_list):
    N = len(x_list)
    delarg = (N*sum(square(x_list)))-(sum(x_list))**2
    A = (sum(square(x_list))*sum(y_list) - sum(x_list)*sum(mult(x_list, y_list)))/delarg
    B = (N*sum(mult(x_list, y_list)) - sum(x_list)*sum(y_list))/delarg
    sig_y = sqrt((1/(N-2))*(sum(square(y_list)) - 2*A*sum(y_list) - 2*B*sum(mult(x_list, y_list)) + 2*A*B*sum(x_list) + B**2*sum(square(x_list)) + N*A**2))
    u_A = sig_y*sqrt(sum(square(x_list))/delarg)
    u_B = sig_y*sqrt(N/delarg)
    print(f'--- LINEAR ---\nLinear (A): {A} +/- {u_A}\nAngular (B): {B} +/- {u_B}')

def ang_coef(x_list, y_list):
    N = len(x_list)
    delarg = (N*sum(square(x_list)))-(sum(x_list))**2
    return (N*sum(mult(x_list, y_list)) - sum(x_list)*sum(y_list))/delarg

def reglin_w(x_list, ux_list, y_list, uy_list):
    ang = ang_coef(x_list, y_list)
    new_uy = []
    for i in range(len(x_list)):
        uy = sqrt(uy_list[i]**2 + ang**2*ux_list[i])
        new_uy.append(uy)
    w_list = inv_square(new_uy)
    sw = sum(w_list)
    swx = sum(mult(w_list, x_list))
    swy = sum(mult(w_list, y_list))
    swxy = sum(mult(mult(w_list, x_list), y_list))
    swx2 = sum(mult(w_list, square(x_list)))
    delarg = sw*swx2 - (swx)**2
    A = (swy*swx2 - swx*swxy)/delarg
    B = (sw*swxy - swx*swy)/delarg
    u_A = sqrt(swx2/delarg)
    u_B = sqrt(sw/delarg)
    print(f'--- PONDERADA ---\nLinear (A): {A} +/- {u_A}\nAngular (B): {B} +/- {u_B}')


reglin_conv(x_list, y_list)
reglin_w(x_list,ux_list,y_list,uy_list)
