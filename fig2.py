from numpy import *
from pylab import *

nt = 129 # of time step
nx = 101 # of grid points
nl = 41 # of layers
xm = 50
tl = xm * 0.75

p  = loadtxt('../base/t.chem.xyz.tsv')

x  = p[0:nx, 0]

tt = p[:, 3]
tr = reshape(tt, (nt, nx*nl))
t  = tr[:, 1]

it1 =  18 # t =  10 day
it2 =  28 # t =  20 day
it3 =  58 # t =  50 day
it4 =  68 # t = 100 day
it5 =  78 # t = 150 day
it6 = 88 # t = 200 day
#it6 = 128 # t = 400 day

delta = 0.5
x = arange(0.0, 50.5, delta)
y = arange(0.0, 20.5, delta)

cc = p[:, 5]
cr = reshape(cc, (nt, nx*nl))
no31 = cr[it1, :]
no32 = cr[it2, :]
no33 = cr[it3, :]
no34 = cr[it4, :]
no35 = cr[it5, :]
no36 = cr[it6, :]

cc = p[:, 10]
cr = reshape(cc, (nt, nx*nl))
u4s1 = cr[it1, :]
u4s2 = cr[it2, :]
u4s3 = cr[it3, :]
u4s4 = cr[it4, :]
u4s5 = cr[it5, :]
u4s6 = cr[it6, :]

for index in range(len(u4s1)):
	u4s1[index] = max(u4s1[index], 1.0001) 
	u4s2[index] = max(u4s2[index], 1.0001) 
	u4s3[index] = max(u4s3[index], 1.0001) 
	u4s4[index] = max(u4s4[index], 1.0001) 
	u4s5[index] = max(u4s5[index], 1.0001) 
	u4s6[index] = max(u4s6[index], 1.0001) 

cc = p[:, 11]
cr = reshape(cc, (nt, nx*nl))
so41 = cr[it1, :]
so42 = cr[it2, :]
so43 = cr[it3, :]
so44 = cr[it4, :]
so45 = cr[it5, :]
so46 = cr[it6, :]

cx = [7.891, 10.000, 10.085, 9.986, 14.429, 19.955, 45.172]
cy = [9.783, 10.000, 8.532, 11.386, 9.244, 8.485, 11.605]
ms = 3

fig, axes = subplots(nrows=3, ncols=5)

subplots_adjust(hspace=0.01, wspace=0.1)

ax = subplot(3, 5, 1)
zm = arange(-1.0e-18, 0.15, 0.01)
z = transpose(reshape(no31, (nl, nx)))
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
title('t = 10')
ylabel('$\mathrm{NO_3^-}$ (mM)')
xticks([])
yticks([])

text(-5, -10, 'Figure 2', fontsize='x-large', color='b')

ax = subplot(3, 5, 2)
z = transpose(reshape(no32, (nl, nx)))
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
title('20')
xticks([])
yticks([])

ax = subplot(3, 5, 3)
z = transpose(reshape(no33, (nl, nx)))
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
title('50')
xticks([])
yticks([])

ax = subplot(3, 5, 4)
z = transpose(reshape(no34, (nl, nx)))
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
title('100')
xticks([])
yticks([])

ax = subplot(3, 5, 5)
z = transpose(reshape(no35, (nl, nx)))
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
title('150 day')
xticks([])
yticks([])

cax1 = fig.add_axes([0.91, 0.642, 0.03, 0.25])
cb=fig.colorbar(c, cax=cax1)

ax = subplot(3, 5, 6)
zmu = arange(0, 4, 0.5)
#zm = arange(-1.0e-20, 1100, 100)
z = transpose(reshape(u4s1, (nl, nx)))
c = contourf(y, x, log10(z), zmu)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
ylabel('U(IV) (log$\mathrm{_{10}}$ $\mathrm{\mu}M$)')
xticks([])
yticks([])

ax = subplot(3, 5, 7)
z = transpose(reshape(u4s2, (nl, nx)))
c = contourf(y, x, log10(z) , zmu)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 8)
z = transpose(reshape(u4s3, (nl, nx)))
c = contourf(y, x, log10(z), zmu)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 9)
z = transpose(reshape(u4s4, (nl, nx)))
c = contourf(y, x, log10(z), zmu)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 10)
z = transpose(reshape(u4s5, (nl, nx)))
c = contourf(y, x, log10(z), zmu)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

cax2 = fig.add_axes([0.91, 0.375, 0.03, 0.25])
cb=fig.colorbar(c, cax=cax2)

ax = subplot(3, 5, 11)
z = transpose(reshape(so41, (nl, nx)))
zm = arange(0, 1.5, 0.1)
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
ylabel('$\mathrm{SO_4^{2-}}$ (mM)')
xticks([])
yticks([])

ax = subplot(3, 5, 12)
z = transpose(reshape(so42, (nl, nx)))
zm = arange(0, 1.5, 0.1)
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 13)
z = transpose(reshape(so43, (nl, nx)))
zm = arange(0, 1.5, 0.1)
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 14)
z = transpose(reshape(so44, (nl, nx)))
zm = arange(0, 1.5, 0.1)
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

ax = subplot(3, 5, 15)
z = transpose(reshape(so45, (nl, nx)))
zm = arange(0, 1.5, 0.1)
c = contourf(y, x, z, zm)
plot(cy, cx, 'ws', markersize=ms, markeredgecolor='w')
axis([0, 20, 50, 0])
ax.set_aspect('equal')
xticks([])
yticks([])

cax3 = fig.add_axes([0.91, 0.108, 0.03, 0.25])
cb=fig.colorbar(c, cax=cax3)

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(8, 11.25)
savefig('fig2.pdf')
show()
