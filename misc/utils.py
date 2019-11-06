# these file have changed location, see /sp4-fwi/data_samples/
rec_core_path = 'output-core/rec.txt'
rec_ops_path  = 'output-ops/rec.txt'

rec_core = []
rec_ops = []

with open(rec_core_path) as rec_core_file, open(rec_ops_path) as rec_ops_file:
    for line_core, line_ops in zip(rec_core_file, rec_ops_file):
            line_core = [float(i) for i in line_core.strip().split(' ')]
            line_ops  = [float(i) for i in line_ops.strip().split(' ')]

            rec_core.append(line_core)
            rec_ops.append(line_ops)

t_rec_core = list(map(list,zip(*rec_core)))
t_rec_ops  = list(map(list,zip(*rec_ops)))

import matplotlib.pyplot as plt
fig1 = plt.figure()

# n-th receptor for both core and ops
CHOOSE_N = 0 # choose which receptor to compare (from 0 to 100)
ax1 = fig1.add_subplot(151)
ax1.plot(t_rec_core[CHOOSE_N], range(len(t_rec_core[0])))
ax1.plot(t_rec_ops[CHOOSE_N], range(len(t_rec_ops[0])))
ax1.xaxis.tick_top()
plt.axis([-5, 5, 0, 716])
plt.ylim(len(t_rec_core[0]), 0)

# difference between Nth receptors (core - ops)
ax2 = fig1.add_subplot(152)
ax2.plot([a - b for a, b in zip(t_rec_core[CHOOSE_N], t_rec_ops[CHOOSE_N])], range(len(t_rec_core[0])), color='red')
ax2.xaxis.tick_top()
plt.yticks([])
plt.axis([-5, 5, 0, 716])

# all receptors core
ax3 = fig1.add_subplot(153)
for i in range(len(t_rec_core)):
    ax3.plot(t_rec_core[i], range(len(t_rec_core[i])))
ax3.xaxis.tick_top()
plt.axis([-5, 5, 0, 716])
plt.yticks([])

# all receptors ops
ax4 = fig1.add_subplot(154)
for i in range(len(t_rec_core)):
    ax4.plot(t_rec_ops[i], range(len(t_rec_ops[i])))
ax4.xaxis.tick_top()
plt.axis([-5, 5, 0, 716])
plt.yticks([])

# all the differences
ax5 = fig1.add_subplot(155)
for i in range(len(t_rec_core)):
    ax5.plot([a - b for a, b in zip(t_rec_core[i], t_rec_ops[i])], range(len(t_rec_core[i])))
ax5.xaxis.tick_top()
plt.axis([-5, 5, 0, 716])
plt.yticks([])

plt.show()


