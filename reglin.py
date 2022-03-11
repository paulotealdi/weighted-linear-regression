from math import sqrt

x_list = [7.8,15.6,23.4,31.2,39.1] #23,42
y_list = [1.2,2.4,3.6,4.7,6] #17,9
ux_list = [0.041,0.041,0.041,0.041,0.041]
uy_list = [0.058,0.058,0.058,0.058,0.058]

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

def ang_coef(x_list, y_list):
    N = len(x_list)
    delarg = (N*sum(square(x_list)))-(sum(x_list))**2
    return (N*sum(mult(x_list, y_list)) - sum(x_list)*sum(y_list))/delarg

def pearson(x_list, y_list):
    x_avg = sum(x_list)/len(x_list)
    y_avg = sum(y_list)/len(y_list)
    num = 0
    den = 0
    x_den = 0
    y_den = 0
    for i in range(len(x_list)):
        num += (x_list[i] - x_avg)*(y_list[i] - y_avg)
    for k in range(len(y_list)):
        y_den += (y_list[k] - y_avg)**2
    for l in range(len(x_list)):
        x_den += (x_list[l] - x_avg)**2
    den = sqrt(x_den*y_den)
    return num/den

def reglin_conv(x_list, y_list):
    N = len(x_list)
    delarg = (N*sum(square(x_list)))-(sum(x_list))**2
    A = (sum(square(x_list))*sum(y_list) - sum(x_list)*sum(mult(x_list, y_list)))/delarg
    B = (N*sum(mult(x_list, y_list)) - sum(x_list)*sum(y_list))/delarg
    sig_y = sqrt((1/(N-2))*(sum(square(y_list)) - 2*A*sum(y_list) - 2*B*sum(mult(x_list, y_list)) + 2*A*B*sum(x_list) + B**2*sum(square(x_list)) + N*A**2))
    u_A = sig_y*sqrt(sum(square(x_list))/delarg)
    u_B = sig_y*sqrt(N/delarg)
    print(f'--- LINEAR ---\nLinear (A): {A} +/- {u_A}\nAngular (B): {B} +/- {u_B}')

def unc_transfer(old_uy, old_ux, ang):
    new_uy = []
    for i in range(len(x_list)):
        uy = sqrt(old_uy[i]**2 + (ang**2)*(old_ux[i]**2))
        new_uy.append(uy)
    return new_uy

def reglin_w(x_list, ux_list, y_list, uy_list, ang=0, k=0):
    if k == 0:
        ang = ang_coef(x_list, y_list)
        
    new_uy = unc_transfer(uy_list, ux_list, ang)
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

    # mesmo processo com o novo coeficiente angular
    ang = B
    new_uy = []
    for i in range(len(x_list)):
        uy = sqrt(uy_list[i]**2 + (ang**2)*(ux_list[i]**2))
        new_uy.append(uy)
    w_list = inv_square(new_uy)
    sw = sum(w_list)
    swx = sum(mult(w_list, x_list))
    swy = sum(mult(w_list, y_list))
    swxy = sum(mult(mult(w_list, x_list), y_list))
    swx2 = sum(mult(w_list, square(x_list)))
    delarg = sw*swx2 - (swx)**2
    nA = (swy*swx2 - swx*swxy)/delarg
    nB = (sw*swxy - swx*swy)/delarg
    nu_A = sqrt(swx2/delarg)
    nu_B = sqrt(sw/delarg)

    # parametro de comparação
    cpar = abs(nA - A)

    if nu_A < cpar:
        # o parametro de comparação precisa ser satisfeito
        print(f'MAIS UMA ITERAÇÃO {k}')
        reglin_w(x_list,ux_list,y_list,uy_list,nB,k+1)

    print(f'--- PONDERADA ---\nLinear (A): {A} +/- {u_A}\nAngular (B): {B} +/- {u_B}')
    print(f'\nIncertezas de y: {new_uy}')


reglin_conv(x_list, y_list)
reglin_w(x_list,ux_list,y_list,uy_list)
print(f'\nPearson = {pearson(x_list, y_list)}')
