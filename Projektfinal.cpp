#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream> 
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <algorithm> 

using namespace std;

struct Vector2D {
    double x, y;

    Vector2D operator + (const Vector2D& other) const {
        return {x + other.x, y + other.y};
    }
    Vector2D operator - (const Vector2D& other) const {
        return {x - other.x, y - other.y};
    }
    Vector2D operator * (double skalar) const {
        return {x * skalar, y * skalar};
    }
    Vector2D operator / (double skalar) const {
        return {x / skalar, y / skalar};
    }
    double length() const {
        return sqrt(x * x + y * y);
    }
    Vector2D normalized() const {
        double len = length();
        return {x / len, y / len};
    }
    double dot(const Vector2D& other) const {
        return x * other.x + y * other.y;
    }
    double distance(const Vector2D& other) const {
        return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }
};

class Particle {
public:
    Vector2D position;
    Vector2D velocity;
    double mass;
    double radius;

    Particle(double x, double y, double vx, double vy, double m, double r)
        : position({x, y}), velocity({vx, vy}), mass(m), radius(r) {}

    void move(double dt) {
        position = position + velocity * dt;
    }

    void collidewall(double boxsize) {
        if (position.x - radius < 0 || position.x + radius > boxsize) {
            velocity.x = -velocity.x; 
            position.x = max(radius, min(position.x, boxsize - radius));
        }

        if (position.y - radius < 0 || position.y + radius > boxsize) {
            velocity.y = -velocity.y; 
            position.y = max(radius, min(position.y, boxsize - radius));
        }
    }

    double kineticEnergy() const {
        return 0.5 * mass * velocity.length() * velocity.length();
    }
    
    double speed() const {
        return velocity.length();
    }
};

void ResolveCollision(Particle& p1, Particle& p2) {
    Vector2D deltadist = p1.position - p2.position;
    double dist = deltadist.length();

    if (dist < p1.radius + p2.radius) {
        Vector2D normaldist = deltadist.normalized();
        Vector2D relativeVelocity = p1.velocity - p2.velocity;

        double velocitydotNormal = relativeVelocity.dot(normaldist);

        if (velocitydotNormal  > 0) return;

        double e = 1.0; 
        double j = -(1 + e) * velocitydotNormal  / (1 / p1.mass + 1 / p2.mass);

        Vector2D impulse = normaldist * j;
        p1.velocity = p1.velocity + impulse * (1 / p1.mass);
        p2.velocity = p2.velocity - impulse * (1 / p2.mass);

        double overlap = 0.5 * (p1.radius + p2.radius - dist);
        p1.position = p1.position + normaldist * overlap;
        p2.position = p2.position - normaldist * overlap;
    }
}

class Simulation {
public:
    vector<Particle> particles;
    double boxsize;
    double timestep;
    ofstream outputFile; 

    Simulation(int particlesAmount, double size, double dt)
        : boxsize(size), timestep(dt) {
        initializeParticles(particlesAmount);
        outputFile.open("simulation_data.csv"); 
        outputFile << "Time,TotalKineticEnergy,TotalMomentum,Speeds\n"; 
    }

    ~Simulation() {
        if (outputFile.is_open()) {
            outputFile.close(); 
        }
    }

    void initializeParticles(int particlesAmount) {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> distPos(0, boxsize);
        uniform_real_distribution<> distVel(-1, 1);
        uniform_real_distribution<> distMass(1, 5);
        uniform_real_distribution<> distRadius(1, 3);

        for (int i = 0; i < particlesAmount; ++i) {
            double x = distPos(gen);
            double y = distPos(gen);
            double vx = distVel(gen);
            double vy = distVel(gen);
            double mass = distMass(gen);
            double radius = distRadius(gen);
            particles.emplace_back(x, y, vx, vy, mass, radius);
        }
    }

    void sweepandprune() {
        vector<pair<double, int>> particlesX;
        for (int i = 0; i < particles.size(); ++i) {
            particlesX.push_back({particles[i].position.x, i});
        }

        sort(particlesX.begin(), particlesX.end());

        for (int i = 0; i < particlesX.size(); ++i) {
            for (int j = i + 1; j < particlesX.size(); ++j) {
                if (particlesX[j].first - particlesX[i].first > (particles[particlesX[i].second].radius + particles[particlesX[j].second].radius)) {
                    break;
                }

                Particle& p1 = particles[particlesX[i].second];
                Particle& p2 = particles[particlesX[j].second];
                if (fabs(p1.position.y - p2.position.y) <= (p1.radius + p2.radius)) {
                    ResolveCollision(p1, p2);
                }
            }
        }
    }

    void update(double currentTime) {
        for (auto& p : particles) {
            p.move(timestep);
            p.collidewall(boxsize);
        }

        sweepandprune();

        if (outputFile.is_open()) {
            double totalEnergy = totalKineticEnergy();
            double totalMomentumMagnitude = totalMomentum().length(); 
            outputFile << currentTime << "," << totalEnergy << "," << totalMomentumMagnitude << ",";
            for (const auto& p : particles) {
                outputFile << p.speed() << " ";
            }
            outputFile << "\n";
        }
    }

    double totalKineticEnergy() const {
        double totalEnergy = 0.0;
        for (const auto& p : particles) {
            totalEnergy += p.kineticEnergy();
        }
        return totalEnergy;
    }

    Vector2D totalMomentum() const {
        Vector2D totalMomentum = {0.0, 0.0};
        for (const auto& p : particles) {
            totalMomentum = totalMomentum + p.velocity * p.mass;
        }
        return totalMomentum;
    }
};

void RenderParticles(const vector<Particle>& particles) {
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    for (const auto& p : particles) {
        draw_list->AddCircleFilled(ImVec2(p.position.x, p.position.y), (float)p.radius, IM_COL32(255, 165, 0, 255));
    }
}

int main() {
    if (!glfwInit()) {
        cerr << "Nie udało się otworzyć glfw" << endl;
        return -1;
    }

    int particlesAmount = 1000;  
    double boxsize = 500.0;  
    double timestep = 0.1;  

    GLFWwindow* window = glfwCreateWindow(800, 800, "Wizualizacja gazu 2D", NULL, NULL);
    if (window == NULL) {
        cerr << "Nie udało się utworzyć okna." << endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    Simulation sim(particlesAmount, boxsize, timestep);
    double currentTime = 0.0; 

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        sim.update(currentTime);
        currentTime += timestep;

        ImGui::Begin("Symulacja Gazu 2D");
        ImGui::SetNextWindowSize(ImVec2(boxsize, boxsize), ImGuiCond_Always);
        RenderParticles(sim.particles);
        ImGui::End();

        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
