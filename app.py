from math import sqrt

def getVelocity(a, e):
    global G, Mj

    return sqrt(G*Mj/a * (1-e)/(1+e))

G           = 6.67e-11
Mj          = 1.898e27                      # Jupiter's mass
Mi          = 8.9319e22                     # Io's mass
Mc          = 1.0759e23                     # Callisto's mass
Me          = 4.7998e22                     # Europa's mass
Mg          = 1.4819e23                     # Ganymede's mass
AU          = 1.5e11
half_hour   = 60*30

i_ap_v      = getVelocity(0.002818 * AU, 0.0041)  # Io's velocity at aphelion
c_ap_v      = getVelocity(0.012585 * AU, 0.0074) # Callisto's velocity at aphelion
e_ap_v      = getVelocity(0.004485 * AU, 0.009)   # Europa's velocity at aphelion
g_ap_v      = getVelocity(0.007155 * AU, 0.0013) # Ganymede's velocity at aphelion

gravconst_i = G*Mi*Mj
gravconst_c = G*Mc*Mj
gravconst_e = G*Me*Mj
gravconst_g = G*Mg*Mj

# Starting conditions
# Io
xi,yi,zi    = 0.002830*AU,0,0
xvi,yvi,zvi = 0,i_ap_v,0

# Callisto
xc,yc,zc    = 0.012680*AU,0,0
xvc,yvc,zvc = 0,c_ap_v,0

# Europa
xe,ye,ze    = 0.004525*AU,0,0
xve,yve,zve = 0,e_ap_v,0

# Ganymede
xg,yg,zg    = 0.007163*AU,0,0
xvg,yvg,zvg = 0,g_ap_v,0

# Jupiter aka sun
xj, yj, zj    = 0,0,0
xvj, yvj, zvj = 0,0,0

t           = 0.0
dt          = 1*half_hour # every frame move this time

xilist,yilist,zilist = [],[],[]
xclist,yclist,zclist = [],[],[]
xelist,yelist,zelist = [],[],[]
xglist,yglist,zglist = [],[],[]
xjlist,yjlist,zjlist = [],[],[]



# start simulation
while t < 1000*365*half_hour:
    ################ Io #############
    # compute G force on Io
    rx, ry, rz = xi - xj, yi - yj, zi - zj
    modr3_i = (rx**2+ry**2+rz**2)**1.5
    fx_i = -gravconst_i*rx/modr3_i
    fy_i = -gravconst_i*ry/modr3_i
    fz_i = -gravconst_i*rz/modr3_i
    
    # update quantities
    xvi += fx_i*dt/Mi
    yvi += fy_i*dt/Mi
    zvi += fz_i*dt/Mi
    
    # update position
    xi += xvi*dt
    yi += yvi*dt
    zi += zvi*dt
    
    # save the position in list
    xilist.append(xi)
    yilist.append(yi)
    zilist.append(zi)
    
    ################ Callistro ##############
    # compute G force on Callistro
    rx_c,ry_c,rz_c = xc - xj, yc - yj, zc - zj
    modr3_c = (rx_c**2+ry_c**2+rz_c**2)**1.5
    fx_c = -gravconst_c*rx_c/modr3_c
    fy_c = -gravconst_c*ry_c/modr3_c
    fz_c = -gravconst_c*rz_c/modr3_c
    
    xvc += fx_c*dt/Mc
    yvc += fy_c*dt/Mc
    zvc += fz_c*dt/Mc
    
    # update position
    xc += xvc*dt
    yc += yvc*dt
    zc += zvc*dt
    
    # add to list
    xclist.append(xc)
    yclist.append(yc)
    zclist.append(zc)
 
    ################ Europa ##############
    # compute G force on Europa
    rx_e,ry_e,rz_e = xe - xj, ye - yj, ze - zj
    modr3_e = (rx_e**2+ry_e**2+rz_e**2)**1.5
    fx_e = -gravconst_e*rx_e/modr3_e
    fy_e = -gravconst_e*ry_e/modr3_e
    fz_e = -gravconst_e*rz_e/modr3_e
    
    xve += fx_e*dt/Me
    yve += fy_e*dt/Me
    zve += fz_e*dt/Me
    
    # update position
    xe += xve*dt
    ye += yve*dt 
    ze += zve*dt
    
    # add to list
    xelist.append(xe)
    yelist.append(ye)
    zelist.append(ze)

    ################ Ganymede ##############
    # compute G force on Ganymede
    rx_g,ry_g,rz_g = xg - xj, yg - yj, zg - zj
    modr3_g = (rx_g**2+ry_g**2+rz_g**2)**1.5
    fx_g = -gravconst_g*rx_g/modr3_g
    fy_g = -gravconst_g*ry_g/modr3_g
    fz_g = -gravconst_g*rz_g/modr3_g
    
    xvg += fx_g*dt/Mg
    yvg += fy_g*dt/Mg
    zvg += fz_g*dt/Mg
    
    # update position
    xg += xvg*dt
    yg += yvg*dt 
    zg += zvg*dt
    
    # add to list
    xglist.append(xg)
    yglist.append(yg)
    zglist.append(zg)
    
    ################ Jupiter ###########
    # update quantities
    xvj += -(fx_i + fx_c + fx_e + fx_g) * dt/Mj
    yvj += -(fx_i + fx_c + fx_e + fx_g) * dt/Mj
    zvj += -(fx_i + fx_c + fx_e + fx_g) * dt/Mj
    
    # update position
    xj += xvj*dt
    yj += yvj*dt
    zj += zvj*dt
    xjlist.append(xj)
    yjlist.append(yj)
    zjlist.append(zj)
    
    # update dt
    t += dt

print('data ready')



import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib

matplotlib.rcParams['animation.embed_limit'] = 2**128

fig, ax = plt.subplots(figsize=(10,10))
ax.set_aspect('equal')
ax.grid(markersize = 0.001)

line_i,     = ax.plot([],[],'-g',lw=1)
point_i,    = ax.plot([0.002830*AU], [0], marker="o", markersize=3, markeredgecolor="orange", markerfacecolor="orange")
text_i      = ax.text(0.002830*AU,0,'Io')

line_c,     = ax.plot([],[],'-g',lw=1)
point_c,    = ax.plot([0.012680*AU], [0], marker="o", markersize=3, markeredgecolor="grey", markerfacecolor="grey")
text_c      = ax.text(0.012680*AU,0,'Callisto')

line_e,     = ax.plot([],[],'-g',lw=1)
point_e,    = ax.plot([0.004525*AU], [0], marker="o", markersize=3, markeredgecolor="blue", markerfacecolor="blue")
text_e      = ax.text(0.004525*AU,0,'Europa')

line_g,     = ax.plot([],[],'-g',lw=1)
point_g,    = ax.plot([0.007163*AU], [0], marker="o", markersize=3, markeredgecolor="black", markerfacecolor="black")
text_g      = ax.text(0.007163*AU,0,'Ganymede')

point_j,    = ax.plot([0], [0], marker="o", markersize=7, markeredgecolor="brown", markerfacecolor="brown")
text_j      = ax.text(0,0,'Jupiter')

ixdata,iydata = [],[]                   # Io track
cxdata,cydata = [],[]                   # Callisto track
exdata,eydata = [],[]                   # Europa track
gxdata,gydata = [],[]                   # Ganymede track
jxdata,jydata = [],[]                   # jupiter track

print(len(xelist))



def update(i):
    ixdata.append(xilist[i])
    iydata.append(yilist[i])
    
    cxdata.append(xclist[i])
    cydata.append(yclist[i])
    
    exdata.append(xelist[i])
    eydata.append(yelist[i])
    
    gxdata.append(xglist[i])
    gydata.append(yglist[i])

    jxdata.append(xjlist[i])
    jydata.append(yjlist[i])
    
    line_i.set_data(ixdata,iydata) 
    point_i.set_data(xilist[i],yilist[i])
    text_i.set_position((xilist[i],yilist[i]))
    
    line_c.set_data(cxdata,cydata)
    point_c.set_data(xclist[i],yclist[i])
    text_c.set_position((xclist[i],yclist[i]))
    
    line_e.set_data(exdata,eydata)
    point_e.set_data(xelist[i],yelist[i])
    text_e.set_position((xelist[i],yelist[i]))

    line_g.set_data(gxdata,gydata)
    point_g.set_data(xglist[i],yglist[i])
    text_g.set_position((xglist[i],yglist[i]))
    
    
    point_j.set_data(xjlist[i],yjlist[i])
    text_j.set_position((xjlist[i],yjlist[i]))
    
    ax.axis('equal')
    ax.set_xlim(-0.015*AU,0.015*AU)
    ax.set_ylim(-0.015*AU,0.015*AU)
    
    return line_i, point_i, text_i, line_c, point_c, text_c, line_e, point_e, text_e, line_g, point_g, text_g, point_j, text_j

anim = animation.FuncAnimation(fig,func=update,frames=len(xelist),interval=1,blit=True)
plt.show()
