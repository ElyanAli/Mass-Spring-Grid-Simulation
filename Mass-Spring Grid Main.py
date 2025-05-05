from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys, math, random

# Global simulation parameters
grid_rows = 10
grid_cols = 10
spacing = 40.0         # distance between adjacent mass points
mass = 1.0             # mass of each point

# 2 marks, set reasonable parameters
k = 500                # spring constant
damping = 1         # damping coefficient for velocity
gravity = -9.8          # gravitational acceleration (downward)
restitution = 0.7      # restitution for collisions

dt = 0.02             # time step (~60 FPS)
# Window and simulation coordinate parameters
window_width = 600
window_height = 800
ortho_left, ortho_right = -300, 300
ortho_bottom, ortho_top = -400, 400  # vertical range changed


class MassSpringGrid:
    def __init__(self, rows, cols, spacing, mass, k, damping, gravity, restitution):
        self.rows = rows
        self.cols = cols
        self.spacing = spacing
        self.mass = mass
        self.k = k
        self.damping = damping
        self.gravity = gravity
        self.restitution = restitution

        self.vertices = []            # list of [x, y] positions
        self.velocities = []          # list of [vx, vy] velocities
        self.springs = []             # list of (i, j) index pairs for connected vertices
        self.spring_rest_lengths = [] # rest length for each spring

        self.init_grid()

    def init_grid(self):
        # Compute starting position so that grid is centered horizontally.
        start_x = - (self.cols - 1) * self.spacing / 2.0
        # Center the grid vertically in (-400, 400)
        start_y = 0.0

        # Create vertices in row-major order.
        for j in range(self.rows):
            for i in range(self.cols):
                x = start_x + i * self.spacing
                y = start_y + j * self.spacing
                self.vertices.append([x, y])

        # Randomize initial velocities.
        self.velocities = [
            [random.uniform(-50, 50), random.uniform(-50, 50)]
            for _ in self.vertices
        ]

        # Define springs connecting neighboring vertices.
        for j in range(self.rows):
            for i in range(self.cols):
                idx = i + j * self.cols
                # Horizontal connection (to right neighbor)
                if i < self.cols - 1:
                    self.springs.append((idx, idx + 1))
                # Vertical connection (to bottom neighbor)
                if j < self.rows - 1:
                    self.springs.append((idx, idx + self.cols))
                # Diagonal: bottom-right neighbor
                if (i < self.cols - 1) and (j < self.rows - 1):
                    self.springs.append((idx, idx + self.cols + 1))
                # Diagonal: bottom-left neighbor
                if (i > 0) and (j < self.rows - 1):
                    self.springs.append((idx, idx + self.cols - 1))

        # Pre-calculate the rest lengths for each spring.
        for (i, j) in self.springs:
            dx = self.vertices[i][0] - self.vertices[j][0]
            dy = self.vertices[i][1] - self.vertices[j][1]
            length = math.sqrt(dx * dx + dy * dy)
            self.spring_rest_lengths.append(length)

    def update_gravity_and_damping_force(self, forces):
        # Add gravity and damping forces.
        for i, (vx, vy) in enumerate(self.velocities):
            # 1 mark
            # Apply gravity force (F = mg)
            forces[i][1] += self.mass * self.gravity

            # Apply damping force (F = -damping * v)
            forces[i][0] -= self.damping * vx
            forces[i][1] -= self.damping * vy
    
    def update_spring_force(self, forces):
        # Calculate spring forces using Hooke’s law.
        for idx, (i, j) in enumerate(self.springs):
            # 3 mark
            # F = -k(L-L0)d
            #positions of the two mass points
            pos_i = self.vertices[i]
            pos_j = self.vertices[j]

            # direction d
            dx = pos_i[0] - pos_j[0]
            dy = pos_i[1] - pos_j[1]

            lengthL = math.sqrt(dx * dx + dy * dy)

            #avoid division by 0
            if lengthL < 0: 
                continue

            dir_x = dx / lengthL
            dir_y = dy / lengthL

            lengthLo = self.spring_rest_lengths[idx]
            #F = -K(L-Lo)
            force_magnitude = -self.k * (lengthL - lengthLo)

            #F = -K(L-Lo) d
            force_x = force_magnitude * dir_x
            force_y = force_magnitude * dir_y

            # Apply to both ends
            forces[i][0] += force_x
            forces[i][1] += force_y
            forces[j][0] -= force_x
            forces[j][1] -= force_y

    def semi_implicit_euler_method(self, selected_index, forces):
        # Update velocities and positions using Euler integration.
        for i in range(len(self.vertices)):
            # 2 mark
            # Skip the selected vertex (being dragged by mouse)
            if i == selected_index:
                continue
            
            #acceleration a = F/m
            acceleration_x = forces[i][0] / self.mass
            acceleration_y = forces[i][1] / self.mass
            
            #Semi-Implicit Euler: v(t+dt) = v(t) + a(t) * dt
            self.velocities[i][0] += acceleration_x * dt
            self.velocities[i][1] += acceleration_y * dt
            
            # Update position p(t+dt) = p(t) + v(t+dt) * dt (1 mark)
            self.vertices[i][0] += self.velocities[i][0] * dt
            self.vertices[i][1] += self.velocities[i][1] * dt


    def collision_hanlding(self):
        # Boundary collision with restitution.
        for i in range(len(self.vertices)):
            # 2 mark
            vert_x = self.vertices[i][0]
            vert_y = self.vertices[i][1]

            velocity_x = self.velocities[i][0] 
            velocity_y = self.velocities[i][1]

            #left boundary check
            if vert_x < ortho_left:
                self.vertices[i][0] = ortho_left
                self.velocities[i][0] = -velocity_x * self.restitution

            #right boundary check
            elif vert_x > ortho_right:
                self.vertices[i][0] = ortho_right
                self.velocities[i][0] = -velocity_x * self.restitution

            #bottom boundary check
            if vert_y < ortho_bottom:
                self.vertices[i][1] = ortho_bottom
                self.velocities[i][1] = -velocity_y * self.restitution

            #top boundary check
            elif vert_y > ortho_top:
                self.vertices[i][1] = ortho_top
                self.velocities[i][1] = -velocity_y * self.restitution



    def update(self, dt, selected_index=None, mouse_pos=None):
        # Create a force accumulator for each vertex.
        forces = [[0.0, 0.0] for _ in self.vertices]

        self.update_gravity_and_damping_force(forces)
        self.update_spring_force(forces)
        # If a vertex is being dragged, override its physics.
        if selected_index is not None and mouse_pos is not None:
            i = selected_index
            self.vertices[i][0] = mouse_pos[0]
            self.vertices[i][1] = mouse_pos[1]
            self.velocities[i] = [0.0, 0.0]

        self.semi_implicit_euler_method(selected_index, forces)
        self.collision_hanlding()

       
    def draw(self):
        # Draw springs as white lines.
        glColor3f(1.0, 1.0, 1.0)
        glBegin(GL_LINES)
        for (i, j) in self.springs:
            glVertex2f(self.vertices[i][0], self.vertices[i][1])
            glVertex2f(self.vertices[j][0], self.vertices[j][1])
        glEnd()

        # Draw mass points as red dots.
        glPointSize(4)
        glBegin(GL_POINTS)
        glColor3f(1.0, 0.0, 0.0)
        for v in self.vertices:
            glVertex2f(v[0], v[1])
        glEnd()


class Simulation:
    def __init__(self):
        self.window_width = window_width
        self.window_height = window_height
        self.ortho_left = ortho_left
        self.ortho_right = ortho_right
        self.ortho_bottom = ortho_bottom
        self.ortho_top = ortho_top

        self.grid = MassSpringGrid(
            grid_rows, grid_cols, spacing, mass, k, damping, gravity, restitution
        )

        self.selected_index = None  # currently dragged vertex index
        self.mouse_sim_pos = (0, 0)   # mouse position in simulation coordinates

    def window_to_sim(self, x, y):
        """Convert window coordinates (origin at top-left) to simulation coordinates."""
        sim_x = self.ortho_left + (x / self.window_width) * (self.ortho_right - self.ortho_left)
        sim_y = self.ortho_top - (y / self.window_height) * (self.ortho_top - self.ortho_bottom)
        return (sim_x, sim_y)

    def display(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        self.grid.draw()
        glutSwapBuffers()

    def update(self, value):
        self.grid.update(dt, self.selected_index, self.mouse_sim_pos)
        glutPostRedisplay()
        glutTimerFunc(int(dt * 50), self.update, 0)

    def mouse_click(self, button, state, x, y):
        sim_pos = self.window_to_sim(x, y)
        if button == GLUT_LEFT_BUTTON:
            if state == GLUT_DOWN:
                threshold = 10.0
                closest = None
                min_dist = float('inf')
                for idx, v in enumerate(self.grid.vertices):
                    dx = v[0] - sim_pos[0]
                    dy = v[1] - sim_pos[1]
                    dist = math.sqrt(dx * dx + dy * dy)
                    if dist < threshold and dist < min_dist:
                        min_dist = dist
                        closest = idx
                self.selected_index = closest
                if self.selected_index is not None:
                    self.mouse_sim_pos = sim_pos
            elif state == GLUT_UP:
                self.selected_index = None

    def mouse_motion(self, x, y):
        self.mouse_sim_pos = self.window_to_sim(x, y)

    def init_gl(self):
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(self.ortho_left, self.ortho_right,
                self.ortho_bottom, self.ortho_top, -1, 1)
        glMatrixMode(GL_MODELVIEW)

    def run(self):
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(self.window_width, self.window_height)
        glutCreateWindow("Mass-Spring Grid Simulation")
        self.init_gl()
        glutDisplayFunc(self.display)
        glutTimerFunc(int(dt * 50), self.update, 0)
        glutMouseFunc(self.mouse_click)
        glutMotionFunc(self.mouse_motion)
        glutMainLoop()


if __name__ == "__main__":
    sim = Simulation()
    sim.run()
