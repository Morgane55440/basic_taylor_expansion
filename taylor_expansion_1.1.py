from copy import deepcopy
import math

def polynomial_approx(E_list_i,n,a):
    E_list=deepcopy(E_list_i)
    factors=(n+1)*[1]
    for k in range(n+1):
        factors[k]=Compute(E_list,a)/(math.factorial(k))
        E_list=derivate(E_list)
        for i in range(k):
            factors[i]+=(factors[k])*((-a)**(k-i))*((math.factorial(k))/(math.factorial(i)*math.factorial(k-i)))
    polynomial=["sum",[]]
    for k in range(len(factors)):
        polynomial[1].append(["product",[factors[k],["power",["x",k]]]])
    return Simplify(polynomial)
    




"""
what we need:  (derivative)
e_sum(a,b,c,....)    (e_sum(da,db,dc,....) )
e_product(a,b,c,....)    (e_sum(e_product(da,b,c,....),e_product(a,db,c,....),e_product(a,b,dc,....),.....))
e_pow(u,c)"u^c"    (e_product(du,c,e_pow(u,e_sum(c,-1))))
specifics:  (derivative)
    e_exp(u)    (e_product(du,e_exp(u)))
    e_ln(u)    (e_product(du,e_pow(u,-1)))
    e_sin(u)    (e_product(du,e_cos(u)))
    e_cos(u)    (e_product(du,e_product(-1,e_sin(u))))
    tan(u)    
    log(a,u)
    asin(u)
    acos(u)
    atan(u)
    sinh(u)
    cosh(u)
    tanh(u)
    asinh(u)
    acosh(u)
    atanh(u)
"""

"""
[func,[elements]]
example:
    5x^2+3x-7:
        ["sum",[-7,["product",["x",3]],["product",[["e_pow",["x",2]],5]]]]
    sin(2x)*(e**(x^2))+(x-3)*(x+1):
        ["sum",[["product",[["e_sin",["product",["x",2]]],["e_exp",["e_pow",["x",2]]]]],["product",[["sum",["x",-3]],["sum",["x",1]]]]]]
"""




def derivate(E_list_i):
    E_list=deepcopy(E_list_i)
    E_list=Simplify(E_list)
    if E_list=="x":
        return(1)
    if isconstant(E_list):
        return(0)
    du=[]
    if E_list[0]=="sum":
        du=["sum",[0]]  #le zéro serat supprimé par Simplify
        for term in E_list[1]:
            dx=derivate(term)
            if dx!=0:
                du[1].append(dx)
        return Simplify(du)
    if E_list[0]=="product":
        du=["sum",[0]] #le zéro sera supprimé par Simplify
        for k in range(len(E_list[1])):
            dx=derivate(E_list[1][k])
            if dx!=0:
                phantom=deepcopy(E_list[1])
                phantom.pop(k)
                phantom.append(dx)
                du[1].append(["product",phantom])
        return Simplify(du)
    if E_list[0]=="power":
        return Simplify(["product",[["power",[E_list[1][0],["sum",[E_list[1][1],-1]]]],derivate(E_list[1][0]),E_list[1][1]]])
    if E_list[0]=="exp":
        return Simplify(["product",[derivate(E_list[1]),E_list]])
    if E_list[0]=="cos":
        return Simplify(["product",[derivate(E_list[1]),["sin",E_list[1]],-1]])
    if E_list[0]=="sin":
        return Simplify(["product",[derivate(E_list[1]),["cos",E_list[1]]]])
    if E_list[0]=="ln":
        return Simplify(["product",[derivate(E_list[1]),["power",[E_list[1],-1]]]])




def term_simplification(term_list_i):
    term_list=deepcopy(term_list_i)
    for k in range(len(term_list)):
        Help_print("si1 : ",term_list[k])
        term_list[k]=Simplify(term_list[k])
        Help_print("si2 : ",term_list[k])
    return term_list




def Simplify_sum_constant(term_list_i):
    term_list=deepcopy(term_list_i)
    c_sum=0
    l=len(term_list)
    Help_print("ssc : ",term_list)
    for k in range(l):
        if isconstant(term_list[k-l]):
            c_sum+=term_list[k-l]
            term_list.pop(k-l)
    if c_sum!=0:
        term_list.append(c_sum)
    return term_list




def Simplify_sum_dupes(term_list_i):  #prend comme entrée une suite de termes dont seul le dernier est une constante
    term_list=deepcopy(term_list_i)
    analysis=[]
    Help_print("ssd : ",term_list)
    for term in term_list:
        found=False
        for dupe_group in analysis:
            if dupe_group[0]==term:
                dupe_group[1]+=1
                found=True
        if not found:
            analysis.append([term,1])
    Help_print("sda : ",analysis)
    for dupe_group in analysis:
        if dupe_group[1]>1:
            for k in range(dupe_group[1]):
                term_list.remove(dupe_group[0])
            term_list.append(["product",dupe_group])
    Help_print("sde : ",term_list)
    return term_list




def constant_factorisation(term_list_i):
    term_list=deepcopy(term_list_i)
    start=deepcopy(term_list)
    Help_print("cfs : ",start)
    done,k=False,0
    l=len(term_list)
    for j in range(l):
        if type(term_list[j])!=list or (type(term_list[j])==list and term_list[j][0]!="product"):
            term_list[j]=["product",[term_list[j],1]]
    while (not done) and (k<l):
        Help_print("cft : ",term_list[k],(type(term_list[k])==list and term_list[k][0]=="product" and isconstant(term_list[k][1][-1])))
        if type(term_list[k])==list and term_list[k][0]=="product" and isconstant(term_list[k][1][-1]):   #potentiellement redundant
            for i in range(k+1,l):
                if type(term_list[i-l])==list and term_list[i-l][0]=="product" and isconstant(term_list[i-l][1][-1]) and unordered_eq(term_list[i-l][1][:-1],term_list[k][1][:-1]):
                    term_list[k][1][-1]+=term_list[i-l][1][-1]
                    term_list.pop(i-l)
                    done=True
        k+=1
    l=len(term_list)
    for j in range(l):
        if type(term_list[j-l])==list and term_list[j-l][0]=="product":
            if len(term_list[j-l][1])==2 and term_list[j-l][1][1]==1:
                term_list[j-l]=term_list[j-l][1][0]
            elif term_list[j-l][1][-1]==0:   # le elif est dût au dernier [1] qui ressort out of range si la première acion a été effectuée et term_list[j-l][1][0] est un seul élément.
                term_list.pop(j)
    Help_print("cfa : ",term_list)
    if start==term_list:
        return start
    return constant_factorisation(term_list)




def Simplify_product_constant(factor_list_i):
    factor_list=deepcopy(factor_list_i)
    Help_print("spc : ",factor_list)
    c_product=1
    l=len(factor_list)
    for k in range(l):
        if isconstant(factor_list[k-l]):
            c_product*=factor_list[k-l]
            factor_list.pop(k-l)
    Help_print("sp2 : ",c_product)
    if c_product!=0:
        if c_product!=1 or len(factor_list)<1:
            factor_list.append(c_product)
        return factor_list
    return []




def Simplify_product_dupes(factor_list_i):
    factor_list=deepcopy(factor_list_i)
    analysis=[]
    Help_print("spd : ",factor_list)
    for factor in factor_list:
        found=False
        for dupe_group in analysis:
            if dupe_group[0]==factor:
                dupe_group[1]+=1
                found=True
        if not found:
            analysis.append([factor,1])
    Help_print("sda : ",analysis)
    for dupe_group in analysis:
        if dupe_group[1]>1:
            for k in range(dupe_group[1]):
                factor_list.remove(dupe_group[0])
            factor_list.append(["power",dupe_group])
    Help_print("sde : ",factor_list)
    return factor_list




def base_grouping(factor_list_i):
    factor_list=deepcopy(factor_list_i)
    start=deepcopy(factor_list)
    Help_print("bss : ",start)
    done,k=False,0
    l=len(factor_list)
    for j in range(l):
        if type(factor_list[j])!=list or (type(factor_list[j])==list and factor_list[j][0]!="power"):
            factor_list[j]=["power",[factor_list[j],1]]
            Help_print("fli : ",factor_list[j])
    Help_print("tep : ",factor_list)
    while (not done) and (k<l):
        if type(factor_list[k])==list and factor_list[k][0]=="power" and isconstant(factor_list[k][1][1]):  #potentiellement redundant
            for i in range(k+1,l):
                if type(factor_list[i-l])==list and factor_list[i-l][0]=="power" and isconstant(factor_list[i-l][1][1]) and unordered_eq(factor_list[i-l][1][0],factor_list[k][1][0]):
                    factor_list[k][1][1]+=factor_list[i-l][1][1]
                    factor_list.pop(i-l)
            done=True
        Help_print("tep : ",factor_list)
        k+=1
    Help_print("tep : ",factor_list)
    l=len(factor_list)
    for j in range(l):
        if type(factor_list[j-l])==list and factor_list[j-l][0]=="power":
            if factor_list[j-l][1][1]==1:
                factor_list[j-l]=factor_list[j-l][1][0]
            elif factor_list[j-l][1][1]==0:  #voir le elif dans constant_factorisation
                if len(factor_list)>1:
                    factor_list.pop(j)
                else :
                    factor_list=[1]
    Help_print("cfa : ",factor_list)
    if start==factor_list:
        return start
    return constant_factorisation(factor_list)




def sub_deconstruction(term_list_i,operation_type):
    term_list=deepcopy(term_list_i)
    l=len(term_list)
    sub_terms=[]
    for k in range(l):
        if type(term_list[k-l])==list and term_list[k-l][0]==operation_type:
            sub_terms=sub_terms+term_list.pop(k-l)[1]
    Help_print("tls : ",term_list,sub_terms)
    term_list+=sub_terms
    if sub_terms==[]:
        return term_list
    return sub_deconstruction(term_list,operation_type)




def Simplify_exp(factor_list_i):
    factor_list=deepcopy(factor_list_i)
    l=len(factor_list)
    exp_list=[]
    for k in range(l):
        if type(factor_list[k-l])==list and factor_list[k-l][0]=="exp":
            exp_list.append(factor_list[k-l][1])
            factor_list.pop(k-l)
    if len(exp_list)>0:
        if len(exp_list)>1:
            exp_list=[Simplify(["sum",exp_list])]
        factor_list.insert(0,["exp",exp_list[0]]) #le but ici est de garder la constante à la fin du produit, bien qu'il me semble que ce soit innnutile, ça ne coûte rien
    return factor_list




def Simplify(E_list_i):
    E_list=deepcopy(E_list_i)
    if isconstant(E_list) or type(E_list)==str:
        return E_list
    if E_list[0]=="sum":
        term_list=E_list[1]
        term_list=constant_factorisation(Simplify_sum_dupes(Simplify_sum_constant(sub_deconstruction(term_simplification(term_list),"sum"))))
        Help_print("tl1 : ",term_list)
        if len(term_list)>1:
            return["sum",term_list]
        if len(term_list)==1:
            return term_list[0]
        return 0
    if E_list[0]=="product":
        factor_list=E_list[1]
        Help_print("tl2 : ",factor_list)
        factor_list=Simplify_exp(base_grouping(Simplify_product_dupes(Simplify_product_constant(sub_deconstruction(term_simplification(factor_list),"product")))))
        Help_print("tl3 : ",factor_list)
        if len(factor_list)>1:
            return["product",factor_list]
        if len(factor_list)==1:
            return factor_list[0]
        return 0   
    if E_list[0]=="power":
        base,exponent=Simplify(E_list[1][0]),Simplify(E_list[1][1])
        if exponent==0:
            return 1
        if exponent==1:
            return base
        if isconstant(base) and isconstant(exponent):
            return base**exponent
        return ["power",[base,exponent]]
    if E_list[0] in ["exp","ln","sin","cos"]:
        if isconstant(E_list[1]):
            for function in [["exp",math.exp],["ln",math.log],["sin",math.sin],["cos",math.cos]]:
                if E_list[0]==function[0]:
                    return function[1](E_list[1])
        return[E_list[0],Simplify(E_list[1])]




def isconstant(ai):
    return (type(a)==int or type(a)==float)




def unordered_eq(element_1_i,element_2_i):
    element_1=deepcopy(element_1_i)
    element_2=deepcopy(element_2_i)
    Help_print("ueq : ",element_1,element_2)
    if element_1==element_2:
        return True
    if isconstant(element_1) and isconstant(element_2) and float(element_1)==float(element_2):
        return True
    if type(element_1)!=type(element_2):
        return False
    for element in element_1:
        if element in element_2:
            element_2.remove(element)
        else:
            return False
    return (element_2==[])




def text_equation(E_list_i):
    E_list=deepcopy(E_list_i)
    if type(E_list)==str:
        return E_list
    if isconstant(E_list):
        if E_list<0:
            E_list="("+str(E_list)+")"
        return str(E_list)
    if E_list[0]=="sum":
        text="("
        term_list=E_list[1]
        for k in range(len(term_list)):
            text+=text_equation(term_list[k])
            if k+1==len(term_list):
                return text+")"
            text+="+"
    if E_list[0]=="product":
        text="("
        factor_list=E_list[1]
        for k in range(len(factor_list)):
            text+=text_equation(factor_list[k])
            if k+1==len(factor_list):
                return text+")"
            text+="*"
    if E_list[0]=="power":
        return "("+text_equation(E_list[1][0])+"^"+text_equation(E_list[1][1])+")"
    if E_list[0] in ["exp","ln","sin","cos"]:
        return E_list[0]+"("+text_equation(E_list[1])+")"




def Help_print(*text): #best function ever
    if Debugging:
        print(*text)




def Compute(E_list_i,x):
    E_list=deepcopy(E_list_i)
    if E_list=="x":
        return x
    if isconstant(E_list):
        return E_list
    if E_list[0]=="sum":
        terms=[]
        for term in E_list[1]:
            terms.append(Compute(term,x))
        return math.fsum(terms)    #is supposed to "Avoids loss of precision"
    if E_list[0]=="product":
        product=1
        for factor in E_list[1]:
            product*=Compute(factor,x)
        return product
    if E_list[0]=="power":
        return (Compute(E_list[1][0],x))**E_list[1][1]
    for function in [["exp",math.exp],["ln",math.log],["sin",math.sin],["cos",math.cos]]:
        if E_list[0]==function[0]:
            return function[1](Compute(E_list[1],x))




def Develop_polynomial(polynomial_i):
    polynomial=deepcopy(polynomial)
    if type(polynomial)==list and polynomial[0]=="sum":
        for k in range(len(polynomial[1])):
            if polynomial[1][k][0]=="product":
                a
                
                #WTF is happening here
        
                
                
                                        
                
                    


    
###########################################################################################################################
Debugging=False
"""        
                
for equation in ["x",5,["sum",[5,"x",-5]] , ["sum",[5,"x",-12,"x"]],["product",[5,"x",["sum",[2,-2]]]] ,
                ["sum",[["sum",["x",["product",["x",2]]]],2,-5,["product",["x",2,0.5,"x"]],["product",[-1,3,"x"]]]],
                ["product",[["sum",["x",4]],["sum",[-5,"x"]],6,-3]],["sum",[25,["product",[-1,["sum",[-12,["product",[7,["sum",[2,["product",[5,4]]]]]]]]]]]],
                 ["product",[["sum",["x",-5]],["power",["x",3]]]],["exp",["power",["x",2]]],["product",[5,["exp",["sin","x"]],["power",[["cos","x"],-1]]]],
                 ["product",[-1,["sin",0]]],["product",[["cos","x"],["exp",2]]]]:
    a=equation
    u=Simplify(a)
    du=derivate(u)
    print("\noriginal equation : ",a,"\nsimplified equation : ",u,"\nderivative : ",du)
    print("original equation : ",text_equation(a),"\nsimplified equation : ",text_equation(u),"\nderivative : ",text_equation(du))
    print("Image of 0 , 1 and 10 : ")
    print("Base equation/Simplified : ",Compute(a,0),"/",Compute(u,0)," ; ",Compute(a,1),"/",Compute(u,1)," ; ",Compute(a,10),"/",Compute(u,10))
    print("Derivate : ",Compute(du,0)," ; ",Compute(du,1)," ; ",Compute(du,10),"\n")
    print("factorisation : ",text_equation(Develop(u)))
"""
for equation in[["square root",["power",["x",0.5]],12,1],["exponential",["exp","x"],8,0],["5",5,10,4],["sin(x)*exp(x)",["product",[["sin","x"],["exp","x"]]],8,0],
                ["ln(x)*cos(x)*exp(x)",["product",[["cos","x"],["exp","x"],["ln","x"]]],9,2],["ln(x)*cos(x)*x",["product",[["cos","x"],"x",["ln","x"]]],9,3.601]]:
               print("The polynomial approximation of degree",equation[2],"of the function which to x affects",
                     equation[0],"around",equation[3],"is :\n",text_equation(polynomial_approx(equation[1],equation[2],equation[3])))    

"""        
for equation in["x",5,["sum",["x",5]],["sum",[4,7,-3,2]],["sum",["x","x"]],["sum",["x",5,3,"x",-12]],["sum",[["sum",[2,"x"]],"x"]],["sum",[["sum",[2,"x"]],["sum",["x"]]]]]:
    print(derivate(Simplify(equation)))
#must return 1,0,1,0,2,2,2,2  
"""        
