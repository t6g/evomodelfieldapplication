from numpy import *
from pylab import *

def myinterp(xi, yi, dx, dy, ic, nt, nxy, p): 
	ix = int(xi/dx)
	iy = int(yi/dy)
	x1 = ix * dx
	x2 = x1 + dx
	y1 = iy * dy
	y2 = y1 + dy
	r11 = (x2 - xi) * (y2 - yi)/dx/dy
	r21 = (xi - x1) * (y2 - yi)/dx/dy
	r12 = (x2 - xi) * (yi - y1)/dx/dy
	r22 = (xi - x1) * (yi - y1)/dx/dy

	cc = p[:, ic]
	cr = reshape(cc, (nt, nxy))
	q11 = cr[:, iy*nx+ix] 
	q21 = cr[:, iy*nx+ix+1]
	q12 = cr[:, (iy+1) * nx + ix]
	q22 = cr[:, (iy+1) * nx + ix + 1]
	z = q11 * r11 + q21 * r21 + q12 * r12 + q22 * r22 

	return z

nt = 129 # of time step
nx = 101 # of grid points
nl = 41 # of layers
xm = 41
nxy = nx * nl

dx = 0.5
dy = 0.5

#FW212 ( 5.300, 10.000)
#MLSC  ( 9.729,  9.224)
#DP13  (15.255,  8.485)
#no = (nl - 1)* nx/2

p  = loadtxt('../base/t.chem.xyz.tsv')

tt = p[:, 3]
tr = reshape(tt, (nt, nxy))
t  = tr[:, 1]

#FW215
xc = 7.891
yc = 9.783
no3fw215 = myinterp(xc, yc, dx, dy,  5, nt, nxy, p) 
fe2fw215 = myinterp(xc, yc, dx, dy,  6, nt, nxy, p) 
u6afw215 = myinterp(xc, yc, dx, dy,  8, nt, nxy, p) 
so4fw215 = myinterp(xc, yc, dx, dy, 11, nt, nxy, p) 
acefw215 = myinterp(xc, yc, dx, dy, 12, nt, nxy, p) 
evofw215 = myinterp(xc, yc, dx, dy, 22, nt, nxy, p) 
evsfw215 = myinterp(xc, yc, dx, dy, 25, nt, nxy, p) 
nrafw215 = myinterp(xc, yc, dx, dy, 13, nt, nxy, p) 
frafw215 = myinterp(xc, yc, dx, dy, 14, nt, nxy, p) 
srofw215 = myinterp(xc, yc, dx, dy, 15, nt, nxy, p) 
srafw215 = myinterp(xc, yc, dx, dy, 16, nt, nxy, p) 
hggfw215 = myinterp(xc, yc, dx, dy, 17, nt, nxy, p) 
megfw215 = myinterp(xc, yc, dx, dy, 18, nt, nxy, p) 

#FW212
xc = 10.0
yc = 10.0
no3fw212 = myinterp(xc, yc, dx, dy,  5, nt, nxy, p) 
fe2fw212 = myinterp(xc, yc, dx, dy,  6, nt, nxy, p) 
u6afw212 = myinterp(xc, yc, dx, dy,  8, nt, nxy, p) 
so4fw212 = myinterp(xc, yc, dx, dy, 11, nt, nxy, p) 
acefw212 = myinterp(xc, yc, dx, dy, 12, nt, nxy, p) 
evofw212 = myinterp(xc, yc, dx, dy, 22, nt, nxy, p) 
evsfw212 = myinterp(xc, yc, dx, dy, 25, nt, nxy, p) 
nrafw212 = myinterp(xc, yc, dx, dy, 13, nt, nxy, p) 
frafw212 = myinterp(xc, yc, dx, dy, 14, nt, nxy, p) 
srofw212 = myinterp(xc, yc, dx, dy, 15, nt, nxy, p) 
srafw212 = myinterp(xc, yc, dx, dy, 16, nt, nxy, p) 
hggfw212 = myinterp(xc, yc, dx, dy, 17, nt, nxy, p) 
megfw212 = myinterp(xc, yc, dx, dy, 18, nt, nxy, p) 


#MLSC
xc = 14.429
yc = 9.224
no3mlsc = myinterp(xc, yc, dx, dy,  5, nt, nxy, p) 
fe2mlsc = myinterp(xc, yc, dx, dy,  6, nt, nxy, p) 
u6amlsc = myinterp(xc, yc, dx, dy,  8, nt, nxy, p) 
so4mlsc = myinterp(xc, yc, dx, dy, 11, nt, nxy, p) 
acemlsc = myinterp(xc, yc, dx, dy, 12, nt, nxy, p) 
evomlsc = myinterp(xc, yc, dx, dy, 22, nt, nxy, p) 
evsmlsc = myinterp(xc, yc, dx, dy, 25, nt, nxy, p) 
nramlsc = myinterp(xc, yc, dx, dy, 13, nt, nxy, p) 
framlsc = myinterp(xc, yc, dx, dy, 14, nt, nxy, p) 
sromlsc = myinterp(xc, yc, dx, dy, 15, nt, nxy, p) 
sramlsc = myinterp(xc, yc, dx, dy, 16, nt, nxy, p) 
hggmlsc = myinterp(xc, yc, dx, dy, 17, nt, nxy, p) 
megmlsc = myinterp(xc, yc, dx, dy, 18, nt, nxy, p) 

#DP13
xc = 19.955
yc =  8.485
no3dp13 = myinterp(xc, yc, dx, dy,  5, nt, nxy, p) 
fe2dp13 = myinterp(xc, yc, dx, dy,  6, nt, nxy, p) 
u6adp13 = myinterp(xc, yc, dx, dy,  8, nt, nxy, p) 
so4dp13 = myinterp(xc, yc, dx, dy, 11, nt, nxy, p) 
acedp13 = myinterp(xc, yc, dx, dy, 12, nt, nxy, p) 
evodp13 = myinterp(xc, yc, dx, dy, 22, nt, nxy, p) 
evsdp13 = myinterp(xc, yc, dx, dy, 25, nt, nxy, p) 
nradp13 = myinterp(xc, yc, dx, dy, 13, nt, nxy, p) 
fradp13 = myinterp(xc, yc, dx, dy, 14, nt, nxy, p) 
srodp13 = myinterp(xc, yc, dx, dy, 15, nt, nxy, p) 
sradp13 = myinterp(xc, yc, dx, dy, 16, nt, nxy, p) 
hggdp13 = myinterp(xc, yc, dx, dy, 17, nt, nxy, p) 
megdp13 = myinterp(xc, yc, dx, dy, 18, nt, nxy, p) 

#SEEP2
xc = 45.172
yc = 11.605
no3seep2 = myinterp(xc, yc, dx, dy,  5, nt, nxy, p) 
fe2seep2 = myinterp(xc, yc, dx, dy,  6, nt, nxy, p) 
u6aseep2 = myinterp(xc, yc, dx, dy,  8, nt, nxy, p) 
so4seep2 = myinterp(xc, yc, dx, dy, 11, nt, nxy, p) 
aceseep2 = myinterp(xc, yc, dx, dy, 12, nt, nxy, p) 
evoseep2 = myinterp(xc, yc, dx, dy, 22, nt, nxy, p) 
evsseep2 = myinterp(xc, yc, dx, dy, 25, nt, nxy, p) 
nraseep2 = myinterp(xc, yc, dx, dy, 13, nt, nxy, p) 
fraseep2 = myinterp(xc, yc, dx, dy, 14, nt, nxy, p) 
sroseep2 = myinterp(xc, yc, dx, dy, 15, nt, nxy, p) 
sraseep2 = myinterp(xc, yc, dx, dy, 16, nt, nxy, p) 
hggseep2 = myinterp(xc, yc, dx, dy, 17, nt, nxy, p) 
megseep2 = myinterp(xc, yc, dx, dy, 18, nt, nxy, p) 

FW215 = loadtxt('../../data/etal/fw215.txt');
FW212 = loadtxt('../../data/etal/fw212.txt');
FW213 = loadtxt('../../data/etal/fw213.txt');
FW214 = loadtxt('../../data/etal/fw214.txt');
MLSC3 = loadtxt('../../data/etal/mlsc3.txt');
MLSC4 = loadtxt('../../data/etal/mlsc4.txt');
MLSC5 = loadtxt('../../data/etal/mlsc5.txt');
MLSE3 = loadtxt('../../data/etal/mlse3.txt');
MLSE4 = loadtxt('../../data/etal/mlse4.txt');
MLSE5 = loadtxt('../../data/etal/mlse5.txt');
DP13  = loadtxt('../../data/etal/dp13.txt');
Seep2 = loadtxt('../../data/etal/seep2.txt');

subplots_adjust(hspace=0.001, wspace=0.00000001)
ms = 5
xmax = 410.0;
xtik = [0, 100, 200, 300, 400]
#acetate
ymax = 3.0;
index = 1;

ax01=subplot(7, 5, 1);
plot(FW215[:, 0], FW215[:, index], 'bo', t, acefw215, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(a)');
xticks(xtik)
ylabel('Acetate')
yticks([0, 1.5, 3])
text(100, 3.2, 'FW215')


text(-150, 4, 'Figure 1', fontsize='x-large', color='b')

ax02=subplot(7, 5, 2);
plot(FW212[:, 0], FW212[:, index], 'bo', FW213[:, 0], FW213[:, index], 'rx', FW214[:, 0], FW214[:, index], 'g^', t, acefw212, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(b)');
xticks(xtik)
yticks([0, 1.5, 3])
text(100, 3.2, 'FW212/3/4')

ax03=subplot(7, 5, 3);
plot(MLSC3[:, 0], MLSC3[:, index], 'bo', MLSC4[:, 0], MLSC4[:, index], 'rx', MLSC5[:, 0], MLSC5[:, index], 'g^', t, acemlsc, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(c)');
xticks(xtik)
yticks([0, 1.5, 3])
text(100, 3.2, 'MLSC')

ax04=subplot(7, 5, 4);
plot(DP13[:, 0], DP13[:, index], 'bo', t, acedp13, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(d)');
xticks(xtik)
yticks([0, 1.5, 3])
text(100, 3.2, 'DP13')

ax05=subplot(7, 5, 5);
plot(Seep2[:, 0], Seep2[:, index], 'bo', t, aceseep2, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
lg=legend(('Obs.',  'Model'), numpoints=1, loc=4)
lg.draw_frame(False)
lgx=lg.get_texts()
setp(lgx, fontsize='small')
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(e)');
xticks(xtik)
yticks([0, 1.5, 3])
text(100, 3.2, 'SEEP2')

# nitrate 
ymax = 0.8;
index = 3;

ax06=subplot(7, 5, 6);
plot(FW215[:, 0], FW215[:, index], 'bo', t, no3fw215, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(f)');
xticks(xtik)
yticks([0, 0.4])
ylabel('Nitrate')

ax07=subplot(7, 5, 7);
plot(FW212[:, 0], FW212[:, index], 'bo', FW213[:, 0], FW213[:, index], 'rx', FW214[:, 0], FW214[:, index], 'g^', t, no3fw212, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
yticks([0, 0.2, 0.4, 0.6, 0.8])
text(0.8*xmax, 0.8*ymax, '(g)');
xticks(xtik)

ax08=subplot(7, 5, 8);
plot(MLSC3[:, 0], MLSC3[:, index], 'bo', MLSC4[:, 0], MLSC4[:, index], 'r+', MLSC5[:, 0], MLSC5[:, index], 'g^', t, no3mlsc, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(h)');
xticks(xtik)
yticks([0, 0.2, 0.4, 0.6, 0.8])

ax09=subplot(7, 5, 9);
plot(DP13[:, 0], DP13[:, index], 'bo', t, no3dp13, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(i)');
xticks(xtik)
yticks([0, 0.2, 0.4, 0.6, 0.8])

ax10=subplot(7, 5, 10);
plot(Seep2[:, 0], Seep2[:, index], 'bo', t, no3seep2, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(j)');
xticks(xtik)
yticks([0, 0.2, 0.4, 0.6, 0.8])

# Fe 
ymax = 0.2;
index = 7;

ax11=subplot(7, 5, 11);
plot(FW215[:, 0], FW215[:, index], 'bo', t, fe2fw215, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(k)');
xticks(xtik)
yticks([0, 0.1])
ylabel('Fe')

ax12=subplot(7, 5, 12);
plot(FW212[:, 0], FW212[:, index], 'bo', FW213[:, 0], FW213[:, index], 'rx', FW214[:, 0], FW214[:, index], 'g^', t, fe2fw212, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(l)');
xticks(xtik)

ax13=subplot(7, 5, 13);
plot(MLSC3[:, 0], MLSC3[:, index], 'bo', MLSC4[:, 0], MLSC4[:, index], 'r+', MLSC5[:, 0], MLSC5[:, index], 'g^', t, fe2mlsc, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(m)');
xticks(xtik)

ax14=subplot(7, 5, 14);
plot(DP13[:, 0], DP13[:, index], 'bo', t, fe2dp13, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(n)');
xticks(xtik)

ax15=subplot(7, 5, 15);
plot(Seep2[:, 0], Seep2[:, index], 'bo', t, fe2seep2, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(o)');
xticks(xtik)

# U(6) 
ymax = 15;
index = 5;

ax16=subplot(7, 5, 16);
plot(FW215[:, 0], FW215[:, index], 'bo', t, u6afw215, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(p)');
xticks(xtik)
yticks([0, 5, 10])
ylabel('U')

ax17=subplot(7, 5, 17);
plot(FW212[:, 0], FW212[:, index], 'bo', FW213[:, 0], FW213[:, index], 'rx', FW214[:, 0], FW214[:, index], 'g^', t, u6afw212, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(q)');
xticks(xtik)
yticks([0, 5, 10, 15])

ax18=subplot(7, 5, 18);
plot(MLSC3[:, 0], MLSC3[:, index], 'bo', MLSC4[:, 0], MLSC4[:, index], 'r+', MLSC5[:, 0], MLSC5[:, index], 'g^', t, u6amlsc, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(r)');
xticks(xtik)
yticks([0, 5, 10, 15])

ax19=subplot(7, 5, 19);
plot(DP13[:, 0], DP13[:, index], 'bo', t, u6adp13, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(s)');
xticks(xtik)
yticks([0, 5, 10, 15])

ax20=subplot(7, 5, 20);
plot(Seep2[:, 0], Seep2[:, index], 'bo', t, u6aseep2, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(t)');
yticks([0, 5, 10, 15])

# sulfate 
ymax = 2.0;
index = 4;

ax21=subplot(7, 5, 21);
plot(FW215[:, 0], FW215[:, index], 'bo', t, so4fw215, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(u)');
yticks([0, 0.5, 1, 1.5])
ylabel('Sulfate')

ax22=subplot(7, 5, 22);
plot(FW212[:, 0], FW212[:, index], 'bo', FW213[:, 0], FW213[:, index], 'rx', FW214[:, 0], FW214[:, index], 'g^', t, so4fw212, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(v)');

ax23=subplot(7, 5, 23);
plot(MLSC3[:, 0], MLSC3[:, index], 'bo', MLSC4[:, 0], MLSC4[:, index], 'r+', MLSC5[:, 0], MLSC5[:, index], 'g^', t, so4mlsc, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(w)');

ax24=subplot(7, 5, 24);
plot(DP13[:, 0], DP13[:, index], 'bo', t, so4dp13, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(x)');

ax25=subplot(7, 5, 25);
plot(Seep2[:, 0], Seep2[:, index], 'bo', t, so4seep2, 'b-', markersize=ms,markerfacecolor='w', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(y)');

ymax = 40
ax26=subplot(7, 5, 26);
plot(t, evsfw215, 'b-', t, evofw215, 'r-', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(z)');
yticks([0, 20])
ylabel('EVO')

ax27=subplot(7, 5, 27);
plot(t, evsfw212, 'b-', t, evofw212, 'r-', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
yticks([0, 10, 20, 30, 40])
text(0.8*xmax, 0.8*ymax, '(A)');

ax28=subplot(7, 5, 28);
plot(t, evsmlsc, 'b-', t, evomlsc, 'r-', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(B)');
yticks([0, 10, 20, 30, 40])

ax29=subplot(7, 5, 29);
plot(t, evsdp13, 'b-', t, evodp13, 'r-', linewidth=2);
xlim([0, xmax]);
ylim([0, ymax]);
text(0.8*xmax, 0.8*ymax, '(C)');
yticks([0, 10, 20, 30, 40])

ax30=subplot(7, 5, 30);
plot(t, evsseep2, 'b-', t, evoseep2, 'r-', linewidth=2);
l = legend(('EVOs', 'CaLCFA$\mathrm{_2}$'), loc = 2);
l.draw_frame(False);
x = l.get_texts();
setp(x, fontsize='small');
xlim([0, xmax]);
ylim([0, ymax]);
yticks([0, 10, 20, 30, 40])
text(0.8*xmax, 0.8*ymax, '(D)');

ymin = 1e-3
ymin = 1e-3
ymax = 100
ax31=subplot(7, 5, 31);
semilogy(t, nrafw215, '-b', t, frafw215, '-.r', t, srofw215+srafw215, '--c', t, hggfw215, ':g', t, megfw215, '-m', linewidth=2)
l = legend(('NRB', 'FeRB'), loc = 2, labelspacing=0.08);
l.draw_frame(False);
x = l.get_texts();
setp(x, fontsize='small');
ylabel('$\mathrm{C_5H_7O_2N}$', fontsize='small');
xlim([0, xmax]);
ylim([ymin, ymax]);
yticks([0.001, 0.1, 10])
#ax31.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
text(0.8*xmax, 7, '(E)');
xticks([0, 100, 200, 300])
#ax = gca()
for tick in ax31.yaxis.get_major_ticks():
	tick.label.set_fontsize('x-small')

ax32=subplot(7, 5, 32);
semilogy(t, nrafw212, '-b', t, frafw212, '-.r', t, srofw212+srafw212, '--c', t, hggfw212, ':g', t, megfw212, '-m', linewidth=2)
xlim([0, xmax]);
ylim([ymin, ymax]);
text(0.8*xmax, 7, '(F)');
xticks([0, 100, 200, 300])

ax33=subplot(7, 5, 33);
semilogy(t, nramlsc, '-b', t, framlsc, '-.r', t, sromlsc+sramlsc, '--c', t, hggmlsc, ':g', t, megmlsc, '-m', linewidth=2)
xlim([0, xmax]);
ylim([ymin, ymax]);
text(0.8*xmax, 7, '(G)');
xticks([0, 100, 200, 300])
xlabel('Elapsed time (day)')

ax34=subplot(7, 5, 34);
semilogy(t, nradp13, '-b', t, fradp13, '-.r', t, srodp13+sradp13, '--c', t, hggdp13, ':g', t, megdp13, '-m', linewidth=2)
xlim([0, xmax]);
ylim([ymin, ymax]);
text(0.8*xmax, 7, '(H)');
xticks([0, 100, 200, 300])

ax35=subplot(7, 5, 35);
semilogy(t, sroseep2+sraseep2, '--c', t, hggseep2, ':g', t, megseep2, '-m', t, nraseep2, '-b', t, fraseep2, '-.r', linewidth=2)
l = legend(('SRB', 'FMB', 'MeG'), loc = 2, labelspacing=0.08);
l.draw_frame(False);
x = l.get_texts();
setp(x, fontsize='small');
xlim([0, xmax]);
ylim([ymin, ymax]);
text(0.8*xmax, 7, '(I)');
xticks([0, 100, 200, 300, 400])

xticklabels=ax01.get_xticklabels()+ax02.get_xticklabels()+ax03.get_xticklabels()+ax04.get_xticklabels()+ax05.get_xticklabels()+ax06.get_xticklabels()+ax07.get_xticklabels()+ax08.get_xticklabels()+ax09.get_xticklabels()+ax10.get_xticklabels()+ax11.get_xticklabels()+ax12.get_xticklabels()+ax13.get_xticklabels()+ax14.get_xticklabels()+ax15.get_xticklabels()+ax16.get_xticklabels()+ax17.get_xticklabels()+ax18.get_xticklabels()+ax19.get_xticklabels()+ax20.get_xticklabels()+ax21.get_xticklabels()+ax22.get_xticklabels()+ax23.get_xticklabels()+ax24.get_xticklabels()+ax25.get_xticklabels()+ax26.get_xticklabels()+ax27.get_xticklabels()+ax28.get_xticklabels()+ax29.get_xticklabels()+ax30.get_xticklabels() 

yticklabels=ax02.get_yticklabels()+ax03.get_yticklabels()+ax04.get_yticklabels()+ax05.get_yticklabels()+ax07.get_yticklabels()+ax08.get_yticklabels()+ax09.get_yticklabels()+ax10.get_yticklabels()+ax12.get_yticklabels()+ax13.get_yticklabels()+ax14.get_yticklabels()+ax15.get_yticklabels()+ax17.get_yticklabels()+ax18.get_yticklabels()+ax19.get_yticklabels()+ax20.get_yticklabels()+ax22.get_yticklabels()+ax23.get_yticklabels()+ax24.get_yticklabels()+ax25.get_yticklabels()+ax27.get_yticklabels()+ax28.get_yticklabels()+ax29.get_yticklabels()+ax30.get_yticklabels()+ax32.get_yticklabels()+ax33.get_yticklabels()+ax34.get_yticklabels()+ax35.get_yticklabels() 

setp(xticklabels, visible=False)
setp(yticklabels, visible=False)

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(10, 10)
savefig('fig1.pdf')
show()

