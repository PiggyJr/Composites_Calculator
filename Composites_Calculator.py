#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np
import pandas as pd
import sqlite3
import matplotlib as mpl
from tkinter import *
from tkinter import messagebox

#Function which refreshes the entry boxes when you update the total number of layers
def onReturn(*args):
    global layer_num
    layer_num = int(layers.get())
    if layer_num == 1:
        layer_num = 1
    else:
        layer_num = int(layer_num/2)
    
    clear_orient()
    clear_thickness()
    for x in range(layer_num):
        global orient_entry
        if layer_num == 1:
            orient_entry = Entry(frame_orient, width=5)
            orient_entry.grid(row=1, column=x+5, pady=20, padx=5)
            break
            
        orient_entry = Entry(frame_orient, width=5)
        orient_entry.grid(row=1, column=x+5, pady=20, padx=5)
        
    for x in range (layer_num):
        global thickness_entry
        if layer_num ==1:
            thickness_entry = Entry(frame_thickness, width=5)
            thickness_entry.grid(row=1, column=x+5, pady=20, padx=5)
            break
            
        thickness_entry = Entry(frame_thickness, width=5)
        thickness_entry.grid(row=1, column=x+5, pady=20, padx=5)

#Functions that clear the the entry boxes before updating to the correct number of layers
def clear_orient(*args):
    for widget in frame_orient.winfo_children():
        widget.destroy()
        
def clear_thickness(*args):
    for widget in frame_thickness.winfo_children():
        widget.destroy()

#Functions which get the data from the orientation and thickness entry boxes
def get_orient():
    children_widgets_orient = frame_orient.winfo_children()
    print("Layup Configuration (Symmetric, °):", end=' ')
    for child_widget_orient in children_widgets_orient:
        print(child_widget_orient.get(), end=' ')

def get_thickness():
    children_widgets_thickness = frame_thickness.winfo_children()
    print("\nLayer Thicknesses (Symmetric, mm):", end=' ')
    for child_widget_thickness in children_widgets_thickness:
        print(child_widget_thickness.get(), end=' ')

#Function which outputs all calculated values
def output():
    if mat_click.get() == mat_options[0]:
        k = 0
        E_x, E_y, E_s, nu_x, X_t, X_c, Y_t, Y_c, S_c = df.loc[0, 'E_x'], df.loc[0, 'E_y'], df.loc[0, 'E_s'], df.loc[0, 'nu_x'], df.loc[0, 'X_t'], df.loc[0, 'X_c'], df.loc[0, 'Y_t'], df.loc[0, 'Y_c'], df.loc[0, 'S_c']
    elif mat_click.get() == mat_options[1]:
        k = 1
        E_x, E_y, E_s, nu_x, X_t, X_c, Y_t, Y_c, S_c = df.loc[1, 'E_x'], df.loc[1, 'E_y'], df.loc[1, 'E_s'], df.loc[1, 'nu_x'], df.loc[1, 'X_t'], df.loc[1, 'X_c'], df.loc[1, 'Y_t'], df.loc[1, 'Y_c'], df.loc[1, 'S_c']
    elif mat_click.get() == mat_options[2]:
        k = 2
        E_x, E_y, E_s, nu_x, X_t, X_c, Y_t, Y_c, S_c = df.loc[2, 'E_x'], df.loc[2, 'E_y'], df.loc[2, 'E_s'], df.loc[2, 'nu_x'], df.loc[2, 'X_t'], df.loc[2, 'X_c'], df.loc[2, 'Y_t'], df.loc[2, 'Y_c'], df.loc[2, 'S_c']
    elif mat_click.get() == mat_options[3]:
        k = 3
        E_x, E_y, E_s, nu_x, X_t, X_c, Y_t, Y_c, S_c = df.loc[3, 'E_x'], df.loc[3, 'E_y'], df.loc[3, 'E_s'], df.loc[3, 'nu_x'], df.loc[3, 'X_t'], df.loc[3, 'X_c'], df.loc[3, 'Y_t'], df.loc[3, 'Y_c'], df.loc[3, 'S_c']
    elif mat_click.get() == mat_options[4]:
        k = 4
        E_x, E_y, E_s, nu_x, X_t, X_c, Y_t, Y_c, S_c = df.loc[4, 'E_x'], df.loc[4, 'E_y'], df.loc[4, 'E_s'], df.loc[4, 'nu_x'], df.loc[4, 'X_t'], df.loc[4, 'X_c'], df.loc[4, 'Y_t'], df.loc[4, 'Y_c'], df.loc[4, 'S_c']
    
    #Calculate Q - on axis (GPa)
    nu_y = nu_x*(E_y/E_x)
    m = 1/(1-(nu_x*nu_y))
    Q_xx = m*E_x
    Q_xy = m*nu_y*E_x
    Q_yx = m*nu_x*E_y
    Q_yy = m*E_y
    Q_ss = E_s
    Q = np.matrix([[Q_xx, Q_xy, 0], [Q_yx, Q_yy, 0], [0, 0, Q_ss]])
    
    #Calculate S - on axis (1/TPa)
    S_xx = (1/E_x)*1000
    S_xy = (-nu_y/E_y)*1000
    S_yx = (-nu_x/E_x)*1000
    S_yy = (1/E_y)*1000
    S_ss = (1/E_s)*1000
    S = np.matrix([[S_xx, S_xy, 0], [S_yx, S_yy, 0], [0, 0, S_ss]])

    def get_stress_matrix():
        global off_axis_stress, sigma_1, sigma_2, sigma_6
        #Getting our inpu off-axis stress values
        sigma_1 = float(sigma_1_entry.get())
        sigma_2 = float(sigma_2_entry.get())
        sigma_6 = float(sigma_6_entry.get())
        #Put our sigmas into a matrix
        off_axis_stress = (np.matrix([[sigma_1], [sigma_2], [sigma_6]]))
        print('\n off axis stress (MPa)')
        print(off_axis_stress)
    
    #Off-axis stress-strain relationship matricies calculations
    def stress_strain_transforms():
        #Setting up lists to store our off-axis matricies in
        global Q_off_axis_list, S_off_axis_list, theta, off_axis_strain_list, stress_transform_list, strain_transform_list
        Q_off_axis_list = []
        S_off_axis_list = []
        off_axis_strain_list = []
        stress_transform_list = []
        strain_transform_list = []
        
        #Setting a ply counter as p
        p = 0
        
        #Calculate U values for [Q] off-axis
        Q_U1 = (1/8)*(3*Q_xx+3*Q_yy+2*Q_xy+4*Q_ss)
        Q_U2 = (1/2)*(Q_xx-Q_yy)
        Q_U3 = (1/8)*(Q_xx+Q_yy-2*Q_xy-4*Q_ss)
        Q_U4 = (1/8)*(Q_xx+Q_yy+6*Q_xy-4*Q_ss)
        Q_U5 = (1/8)*(Q_xx+Q_yy-2*Q_xy+4*Q_ss)
        
        #Calculate U values for [S] off-axis
        S_U1 = (1/8)*(3*S_xx+3*S_yy+2*S_xy+S_ss)
        S_U2 = (1/2)*(S_xx-S_yy)
        S_U3 = (1/8)*(S_xx+S_yy-2*S_xy-S_ss)
        S_U4 = (1/8)*(S_xx+S_yy+6*S_xy-S_ss)
        S_U5 = (1/2)*(S_xx+S_yy-2*S_xy+S_ss)
        
        children_widgets_orient = frame_orient.winfo_children()
        for child_widget_orient in children_widgets_orient:
            theta = float(child_widget_orient.get())*(np.pi/180)
            p = p+1
            
            #Calculating off-axis modulus [Q]
            #Calculate Q_ij values
            Q_11 = Q_U1 + np.cos(2*theta)*Q_U2 + np.cos(4*theta)*Q_U3
            Q_22 = Q_U1 + (-np.cos(2*theta)*Q_U2) + np.cos(4*theta)*Q_U3
            Q_12 = Q_U4 + (-np.cos(4*theta)*Q_U3)
            Q_66 = Q_U5 + (-np.cos(4*theta)*Q_U3)
            Q_16 = (1/2)*np.sin(2*theta)*Q_U2 + np.sin(4*theta)*Q_U3
            Q_26 = (1/2)*np.sin(2*theta)*Q_U2 + (-np.sin(4*theta)*Q_U3)
            #Creating off-axis matrix
            Q_off_axis = np.matrix([[Q_11, Q_12, Q_16],
                                    [Q_12, Q_22, Q_26],
                                    [Q_16, Q_26, Q_66]])
            #Storing our Q_off_axis matrices for each layer in a list
            Q_off_axis_list.append(Q_off_axis)
            
            #Calculating off-axis complience [S]
            #Calculating S_ij values
            S_11 = S_U1 + np.cos(2*theta)*S_U2 + np.cos(4*theta)*S_U3
            S_22 = S_U1 + (-np.cos(2*theta))*S_U2 + np.cos(4*theta)*S_U3
            S_12 = S_U4 + (-np.cos(4*theta))*S_U3
            S_66 = S_U5 + (-4*np.cos(4*theta))*S_U3
            S_16 = np.sin(2*theta)*S_U2 + 2*np.sin(4*theta)*S_U3
            S_26 = np.sin(2*theta)*S_U2 + (-2*np.sin(4*theta))*S_U3
            #Creating off-axis matrix
            S_off_axis = np.matrix([[S_11, S_12, S_16],
                                    [S_12, S_22, S_26],
                                    [S_16, S_26, S_66]])
            #Storing our S_off_axis matricies for each layer in a list
            S_off_axis_list.append(S_off_axis)
            
            #Off-axis strain calculations
            off_axis_strain = (np.matmul(S_off_axis_list[p-1], off_axis_stress))/1000
            #Storing the off-axis strain
            off_axis_strain_list.append(off_axis_strain)
            
            #On-axis stress and strain transforms
            #Calculating on-axis stresses given off-axis stresses
            sigma_x = (1/2)*(sigma_1+sigma_2)+(1/2)*(sigma_1-sigma_2)*np.cos(2*theta)+sigma_6*np.sin(2*theta)
            sigma_y = (1/2)*(sigma_1+sigma_2)-(1/2)*(sigma_1-sigma_2)*np.cos(2*theta)-sigma_6*np.sin(2*theta)
            sigma_s = (-1/2)*(sigma_1-sigma_2)*np.sin(2*theta)+sigma_6*np.cos(2*theta)
            stress_transform = np.matrix([[sigma_x], [sigma_y], [sigma_s]])
            stress_transform_list.append(stress_transform)
            
            #Getting our epsilon_1,2,6 values from our off-axis strain matrix
            epsilon_1 = off_axis_strain.item(0)
            epsilon_2 = off_axis_strain.item(1)
            epsilon_6 = off_axis_strain.item(2)
            
            #Calculating on-axis strains given off-axis strains
            epsilon_x = (1/2)*(epsilon_1+epsilon_2)+(1/2)*(epsilon_1-epsilon_2)*np.cos(2*theta)+(1/2)*epsilon_6*np.sin(2*theta)
            epsilon_y = (1/2)*(epsilon_1+epsilon_2)-(1/2)*(epsilon_1-epsilon_2)*np.cos(2*theta)-(1/2)*epsilon_6*np.sin(2*theta)
            epsilon_s = (-1)*(epsilon_1-epsilon_2)*np.sin(2*theta)+epsilon_6*np.cos(2*theta)
            strain_transform = np.matrix([[epsilon_x], [epsilon_y], [epsilon_s]])
            strain_transform_list.append(strain_transform)
            
            print('\noff-axis [Q] and off-axis [S] for layer %d:' %p)
            print('\n[Q] off axis (GPa):')
            print(np.round_(Q_off_axis, decimals=4))
            print('\n[S] off axis (1/TPa):')
            print(np.round_(S_off_axis, decimals=4))
            print('\nOff-axis strain(mm/m):')
            print(np.round_(off_axis_strain, decimals=4))
            print('\nOn-axis stress (MPa):')
            print(np.round_(stress_transform, decimals=4))
            print('\nOn-axis strain (mm/m):')
            print(np.round_(strain_transform, decimals=4))
            
    
    print("----------------input parameters----------------")
    print(df.iloc[[k, 5],:])
    print("\nTotal Number of Layers:" + layers.get())
    get_orient()
    get_thickness()
    print("\nCore Thickness(mm):" + core.get())
    get_stress_matrix()
    print("\n----------------output matricies----------------")
    print('[Q] on-axis(GPa)')
    print(np.round_(Q, decimals=4))
    print('\n[S] on-axis(1/TPa)')
    print(np.round_(S, decimals=4))
    stress_strain_transforms()
    

#Material Properties Table into a dataframe
mat_prop = {
    'Material':['T300/N5208', 'E-glass/epoxy', 'Kev49/epoxy', 'AS/H3501', 'AS4/PEEK', 'Units'],
    'E_x':[181,38.6,76,138,134,'GPa'],
    'E_y':[10.3,8.27,5.5,8.96,8.9, 'GPa'],
    'nu_x':[0.28,0.26,0.34,0.3,0.28, ''],
    'E_s':[7.17,4.14,2.3,7.1,5.1, 'GPa'],
    'X_t':[1500,1062,1400,1447,2130, 'MPa'],
    'X_c':[1500,610,235,1447,1100, 'MPa'],
    'Y_t':[40,31,12,51.7,80, 'MPa'],
    'Y_c':[246,118,53,206,200, 'MPa'],
    'S_c':[68,72,34,93,160, 'MPa']
}
df = pd.DataFrame(mat_prop)

##GUI##
root = Tk()

root.title("Composite Material Analysis Tool")
root.geometry("900x600+20+10")

#Setting up frame for material selection
frame_matsel = LabelFrame(root, text="Material Selection", padx=5, pady=5)
frame_matsel.place(x=50, y=10)

#Setting up frame for total number of layers
frame_layers = LabelFrame(root, text="Total Number of Layers\n"+"(Press Enter to Apply)", padx=5, pady=5)
frame_layers.place(x=80, y=80)

#Setting up frame for core thickness
frame_core = LabelFrame(root, text="Core thickness (mm)", padx=5, pady=5)
frame_core.place(x=85, y=150)

#Setting up frame for orientation inputs
frame_orient = LabelFrame(root, text="Layup Inputs (symmetric)", padx=5, pady=5)
frame_orient.place(x=300, y=10)

#Setting up frame for thickness inputs
frame_thickness = LabelFrame(root, text="Thickness Inputs (symmetric, mm)", padx=5, pady=5)
frame_thickness.place(x=300, y=110)

#Setting up frame for off-axis stress inputs
frame_stress = LabelFrame(root, text = 'Off-axis stress Inputs (MPa)', padx=5, pady=5)
frame_stress.place(x=75, y=220)

#Setting up the dropdown material selection dropdown menu
mat_click = StringVar()
mat_options = [
    "graphite/epoxy (T300/N5208)",
    "Fiberglass (E-glass/epoxy)",
    "Kevlar/epoxy (Kev49/epoxy)",
    "graphite/epoxy (AS/H3501)",
    "graphite/ thermoplastic (AS4/PEEK)"
]

mat_click.set(mat_options[0])

mat_select = OptionMenu(frame_matsel, mat_click, *mat_options)
mat_select.pack()

#Setting up total number of layers input and orientation inputs
layers = Entry(frame_layers, width=10)
layers.bind("<Return>", onReturn)
layers.pack()

#Setting up core thickness input
core = Entry(frame_core, width=10)
core.pack()

#Setting up output using a button
disp_output = Button(root, text='Click to Output to Console', width=22, command=output)
disp_output.place(x=400, y=400)

#Setting up stress input entries
sigma_1_label = Label(frame_stress, text='σ_1')
sigma_1_label.grid(row=0, column=0)
sigma_2_label = Label(frame_stress, text='σ_2')
sigma_2_label.grid(row=1, column=0)
sigma_6_label = Label(frame_stress, text='σ_6')
sigma_6_label.grid(row=2, column=0)

sigma_1_entry = Entry(frame_stress, width=10)
sigma_1_entry.insert(END,'0')
sigma_1_entry.grid(row=0, column=1)
sigma_2_entry = Entry(frame_stress, width=10)
sigma_2_entry.insert(END,'0')
sigma_2_entry.grid(row=1, column=1)
sigma_6_entry = Entry(frame_stress, width=10)
sigma_6_entry.insert(END,'0')
sigma_6_entry.grid(row=2, column=1)

root.mainloop()


# In[ ]:





# In[ ]:

