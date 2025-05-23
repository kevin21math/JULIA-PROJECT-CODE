import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

U_num = np.loadtxt("C:\\Users\\kevin\\Downloads\\VSCODE PROGRAMS\\dis_advection_upwind_central_data.csv", delimiter=',')
U_exact = np.loadtxt("C:\\Users\\kevin\\Downloads\\VSCODE PROGRAMS\\dis_exact_solution.csv", delimiter=',')

Nt_plus1, Nx = U_num.shape
x = np.linspace(0, 2, Nx)

fig, ax = plt.subplots()
line_num, = ax.plot([], [], 'b-', label='Upwind-Central')
line_exact, = ax.plot([], [], 'r--', label='Exact')
ax.set_xlim(0, 2)
ax.set_ylim(-0.1, 1.1)
ax.set_xlabel("x")
ax.set_ylabel("u")
ax.set_title("1D Linear Advection: Upwind-Central vs Exact")
ax.legend()

def init():
    line_num.set_data([], [])
    line_exact.set_data([], [])
    return line_num, line_exact

def update(frame):
    line_num.set_data(x, U_num[frame])
    line_exact.set_data(x, U_exact[frame])
    ax.set_title(f"timestep = {frame}")
    return line_num, line_exact

ani = animation.FuncAnimation(
    fig, update, frames=Nt_plus1,
    init_func=init, blit=True, interval=50
)

ani.save("dis_advection_upwind_central.gif", writer='pillow', fps=20)

fig_final, ax_final = plt.subplots()
ax_final.plot(x, U_num[-1], 'b-', label='Upwind-Central')
ax_final.plot(x, U_exact[-1], 'r--', label='Exact')
ax_final.set_xlim(0, 2)
ax_final.set_ylim(-0.1, 1.1)
ax_final.set_xlabel("x")
ax_final.set_ylabel("u")
ax_final.set_title("Final Time Step: Upwind-Central vs Exact")
ax_final.legend()
fig_final.savefig("dis_advection_upwind_central.png", dpi=300)