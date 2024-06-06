import * as THREE from 'three';
import * as dat from 'dat.gui';
import Stats from 'stats.js';


// Parameters for a viscous fluid: smoothingLength: 0.55, graivity: 2, referenceDensity: 2.9, gasConstant: 0.4, viscosity: -4, particlesMass: 0.3, timeStep: 0.2

const simulationParams = {
    
    boxWidth: 8,                 // Width of the simulation box
    boxHeight: 4,                // Height of the simulation box
    boxDepth: 4,                 // Depth of the simulation box (only used in 3D)

    numberOfParticles: 8,        // Number of particles in each dimension
    particlesMetalness: 0.3,     // Metalness of the particles, higher values will make the particles to be more reflective
    particlesRoughness: 0.3,     // Roughness of the particles, higher values will make the particles to be less smooth



    particlesDimensions: 0.15,   // Diameter of the particles
    particleSpacing: 0.30,       // Initial spacing between the particles
    zSimulation: 7,              // Z position of the camera
    
    
    ThirdDimension: true,        // Set to false to simulate in 2D
    enableReflections: false,    // Set to true to enable real-time reflections (works only when InstancedMesh is disabled)
    InstancedMesh: true,         // Set to true to use an InstancedMesh for rendering the particles
    useGridOptimization: true,   // Set to true to use a grid optimization for the SPH calculations
    showEnvironment: false,      // Set to true to show an environment map

    BoxRestitution: 0.9,         // Restitution of the box walls, higher values will make the particles to bounce more on the walls
    particlesMass: 0.3,          // Mass of the particles, higher values will affect interactions and forces between particles, impacting the simulation dynamics
    smoothingLength: 0.3,        // Radius of the smoothing kernel, higher values will check more neighbors particles for calculations
    referenceDensity: 1,         // Target density for the fluid, higher values will make the particles to be more compressed
    gasConstant: 0.5,            // Gas constant for the equation of state, higher values will make the fluid more incompressible
    viscosity: 1,                // Viscosity of the fluid, higher values will make the relative velocity of the particles to be more similar
    gravity: 8,                  // Gravity force applied to the particles
    timeStep: 0.1                // Simulation time step

};


// Class to represent a particle in the simulation
class Particle {
    constructor(position, velocity, index) {
        this.position = position;
        this.velocity = velocity;
        this.acceleration = new THREE.Vector3(0, 0, 0);
        this.mass = simulationParams.particlesMass;
        this.density = 0;
        this.pressure = 0;
        this.forces = new THREE.Vector3(0, 0, 0);
        this.index = index;

        if (!simulationParams.InstancedMesh) {

            // Create CubeCamera for real-time reflections
            const cubeRenderTarget = new THREE.WebGLCubeRenderTarget(128, {
                format: THREE.RGBFormat,
                generateMipmaps: true,
                minFilter: THREE.LinearMipmapLinearFilter,
            });

            this.cubeCamera = new THREE.CubeCamera(0.1, 1000, cubeRenderTarget);

            this.sphereMaterial = new THREE.MeshStandardMaterial({
                color: 0xffff00,
                metalness: simulationParams.particlesMetalness,
                roughness: simulationParams.particlesRoughness,
                envMap: this.cubeCamera.renderTarget.texture,
            });

            this.mesh = new THREE.Mesh(sphereGeometry, this.sphereMaterial);
            scene.add(this.mesh);
        }
    }
}

// Function to create the GUI for the simulation parameters
function createGui() {
    const gui = new dat.GUI();

    //set default width
    gui.width = 300;
    gui.add(simulationParams, 'timeStep', 0.01, 0.5).name('Time Step');
    gui.add(simulationParams, 'smoothingLength', 0.01, 2.0).name('Smoothing Length');
    gui.add(simulationParams, 'gravity', 0, 15).name('Gravity');
    gui.add(simulationParams, 'referenceDensity', 0.01, 10).name('Reference Density');
    gui.add(simulationParams, 'gasConstant', 0.01, 2).name('Gas Constant');
    gui.add(simulationParams, 'viscosity', 0, 10).name('Viscosity');

    // On particle mass change, loop over all the particles and update their mass
    gui.add({ particlesMass: simulationParams.particlesMass }, 'particlesMass', 0.1, 1).name('Mass').onChange((value) => {
        simulationParams.particlesMass = value;
        for (let i = 0; i < particles.length; i++) {
            particles[i].mass = simulationParams.particlesMass;
        }
    });

    // On particle metalness change, loop over all the particles and update their metalness
    gui.add({ particlesMetalness: simulationParams.particlesMetalness }, 'particlesMetalness', 0, 1).name('Metalness').onChange((value) => {
        simulationParams.particlesMetalness = value;
        if (!simulationParams.InstancedMesh){
            for (let i = 0; i < particles.length; i++) {
                particles[i].sphereMaterial.metalness = simulationParams.particlesMetalness;
            }   
        } else {
            instancedMaterial.metalness = simulationParams.particlesMetalness;
        }
    });

    // On particle roughness change, loop over all the particles and update their roughness
    gui.add({ particlesRoughness: simulationParams.particlesRoughness }, 'particlesRoughness', 0, 1).name('Roughness').onChange((value) => {
        simulationParams.particlesRoughness = value;
        if (!simulationParams.InstancedMesh){
            for (let i = 0; i < particles.length; i++) {
                particles[i].sphereMaterial.roughness = simulationParams.particlesRoughness;
            }   
        } else {
            instancedMaterial.roughness = simulationParams.particlesRoughness;
        }
    });

    gui.add({ zSimulation: simulationParams.zSimulation }, 'zSimulation', 0, 20).name('Z Distance').onChange((value) => {
        simulationParams.zSimulation = value;
        camera.position.z = value;
    });

    // Create a section for the simulation box dimensions and restitution parameters
    const boxFolder = gui.addFolder('Simulation Box Parameters');
    boxFolder.add(simulationParams, 'boxWidth', 1, 10).name('Box Width').onChange(updateBoxGeometry);
    boxFolder.add(simulationParams, 'boxHeight', 1, 10).name('Box Height').onChange(updateBoxGeometry);
    boxFolder.add(simulationParams, 'boxDepth', 1, 10).name('Box Depth').onChange(updateBoxGeometry);
    boxFolder.add(simulationParams, 'BoxRestitution', 0, 1).name('Box Restitution');
    boxFolder.open();

    // Toggle environment map
    gui.add({ toggleEnvironment: () => {
        simulationParams.showEnvironment = !simulationParams.showEnvironment;
        if (simulationParams.showEnvironment) {
            loader.load('map.jpg', function (texture) {
                texture.mapping = THREE.EquirectangularReflectionMapping;
                scene.background = texture;
                scene.environment = texture;
            });

            // Set the environment map to the particles
            for (let i = 0; i < particles.length; i++) {
                particles[i].sphereMaterial.envMap = scene.environment;
            }
        } else {
            scene.environment = null;
            scene.background = null;
            for (let i = 0; i < particles.length; i++) {
                particles[i].sphereMaterial.envMap = null;
            }
        }
    }}, 'toggleEnvironment').name('Toggle Environment Map');

    gui.add({ toggleGravity: () => {
        simulationParams.gravity = simulationParams.gravity === 0 ? 9 : 0;
        gui.__controllers.find(c => c.property === 'gravity').updateDisplay();
    }}, 'toggleGravity').name('Toggle Gravity');

    gui.add({ toggleReflections: () => {
        if (!simulationParams.InstancedMesh) {
            simulationParams.enableReflections = !simulationParams.enableReflections;
            computeReflections(particles); // Ensure reflections are computed when toggled
        }
    }}, 'toggleReflections').name('Toggle Reflections (InstancedMesh disabled)');

    gui.add({ toggleGridOptimization: () => {
        simulationParams.useGridOptimization = !simulationParams.useGridOptimization;
    }}, 'toggleGridOptimization').name('Toggle Grid Optimization');
}


// Function to update the geometry of the simulation box according to the parameters selected in the GUI
function updateBoxGeometry() {

    scene.remove(cube);
    if (simulationParams.ThirdDimension) {
        cubeGeometry = new THREE.BoxGeometry(simulationParams.boxWidth, simulationParams.boxHeight, simulationParams.boxDepth);
    } else {
        cubeGeometry = new THREE.BoxGeometry(simulationParams.boxWidth, simulationParams.boxHeight, 0.0);
    }
    cubeMaterial = new THREE.MeshBasicMaterial({ color: 0x000000, wireframe: true });
    cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
    scene.add(cube);
}

// Function to create a grid of particles in the simulation box
function gridOfParticles(x, y, z, boxWidth, boxHeight, boxDepth, spacing) {
    const particles = [];
    let index = 0;
    for (let i = 0; i < x; i++) {
        for (let j = 0; j < y; j++) {
            if (simulationParams.ThirdDimension) {
                for (let k = 0; k < z; k++) {
                    const position = new THREE.Vector3(
                        i * spacing - boxWidth / 2 + Math.random() * 0.01,
                        j * spacing - boxHeight / 2 + Math.random() * 0.01,
                        k * spacing - boxDepth / 2 + Math.random() * 0.01
                    );
                    const velocity = new THREE.Vector3(0, 0, 0);
                    const particle = new Particle(position, velocity, index);
                    particles.push(particle);
                    index++;
                }
            } else {
                const position = new THREE.Vector3(
                    i * spacing - boxWidth / 2 + Math.random() * 0.01,
                    j * spacing - boxHeight / 2 + Math.random() * 0.01,
                    0
                );
                const velocity = new THREE.Vector3(0, 0, 0);
                const particle = new Particle(position, velocity, index);
                particles.push(particle);
                index++;
            }
        }
    }
    return particles;
}

// Constants for the SPH kernel functions to avoid recalculating them every frame
const kernelCoefficients = {
    cubic: 315 / (64 * Math.PI),
    gradient: -45 / Math.PI,
};

// Function to calculate the cubic spline kernel function for given distance and smoothing length
function kernel(r, h) {
    const q = r / h;
    if (q < 1) {
        return kernelCoefficients.cubic * Math.pow(h * h - r * r, 3) / Math.pow(h, 9);
    }
    return 0;
}

// Function to calculate the gradient of the cubic spline kernel function for given distance and smoothing length
function gradientKernel(r, h) {
    const q = r / h;
    if (q < 1) {
        return kernelCoefficients.gradient * Math.pow(h - r, 2) / Math.pow(h, 6);
    }
    return 0;
}

// Function to calculate the density and pressure of the particles
function calculateDensityAndPressure(particles, h, rho0, k) {
    for (let i = 0; i < particles.length; i++) {
        let density = 0;
        for (let j = 0; j < particles.length; j++) {
            const r = particles[i].position.distanceTo(particles[j].position);
            if (r < h) { // Only consider neighbors within smoothing length
                density += particles[j].mass * kernel(r, h);
            }
        }
        particles[i].density = density;
        particles[i].pressure = k * (particles[i].density - rho0);
    }
}

// Function to calculate the forces acting on the particles using a grid optimization
function calculateDensityAndPressureOptimized(particles, h, rho0, k) {
    const cellSize = h;
    const cells = new Map();
    for (let i = 0; i < particles.length; i++) {
        const cellX = Math.floor(particles[i].position.x / cellSize);
        const cellY = Math.floor(particles[i].position.y / cellSize);
        const cellZ = Math.floor(particles[i].position.z / cellSize);
        const cellId = `${cellX}-${cellY}-${cellZ}`;
        if (!cells.has(cellId)) {
            cells.set(cellId, []);
        }
        cells.get(cellId).push(i);
    }
    for (let i = 0; i < particles.length; i++) {
        let density = 0;
        const cellX = Math.floor(particles[i].position.x / cellSize);
        const cellY = Math.floor(particles[i].position.y / cellSize);
        const cellZ = Math.floor(particles[i].position.z / cellSize);
        for (let x = cellX - 1; x <= cellX + 1; x++) {
            for (let y = cellY - 1; y <= cellY + 1; y++) {
                for (let z = cellZ - 1; z <= cellZ + 1; z++) {
                    const cellId = `${x}-${y}-${z}`;
                    if (cells.has(cellId)) {
                        for (const j of cells.get(cellId)) {
                            const r = particles[i].position.distanceTo(particles[j].position);
                            if (r < h) { // Only consider neighbors within smoothing length
                                density += particles[j].mass * kernel(r, h);
                            }
                        }
                    }
                }
            }
        }
        particles[i].density = density;
        particles[i].pressure = k * (particles[i].density - rho0);
    }
}

// Function to calculate the forces acting on the particles using a grid optimization
function calculateForcesOptimized(particles, h, mu, g) {
    const cellSize = h;
    const cells = new Map();
    for (let i = 0; i < particles.length; i++) {
        const cellX = Math.floor(particles[i].position.x / cellSize);
        const cellY = Math.floor(particles[i].position.y / cellSize);
        const cellZ = Math.floor(particles[i].position.z / cellSize);
        const cellId = `${cellX}-${cellY}-${cellZ}`;
        if (!cells.has(cellId)) {
            cells.set(cellId, []);
        }
        cells.get(cellId).push(i);
    }
    for (let i = 0; i < particles.length; i++) {
        let pressureForce = new THREE.Vector3(0, 0, 0);
        let viscosityForce = new THREE.Vector3(0, 0, 0);
        const gravityForce = new THREE.Vector3(0, -g * particles[i].mass, 0);
        particles[i].forces.set(0, 0, 0);
        const cellX = Math.floor(particles[i].position.x / cellSize);
        const cellY = Math.floor(particles[i].position.y / cellSize);
        const cellZ = Math.floor(particles[i].position.z / cellSize);
        for (let x = cellX - 1; x <= cellX + 1; x++) {
            for (let y = cellY - 1; y <= cellY + 1; y++) {
                for (let z = cellZ - 1; z <= cellZ + 1; z++) {
                    const cellId = `${x}-${y}-${z}`;
                    if (cells.has(cellId)) {
                        for (const j of cells.get(cellId)) {
                            if (i === j) continue; // Skip the particle itself
                            const r = particles[i].position.distanceTo(particles[j].position);
                            if (r < h) { // Only consider neighbors within smoothing length
                                const gradient = gradientKernel(r, h);
                                const pressure = (particles[i].pressure + particles[j].pressure) / (2 * particles[j].density);
                                const direction = new THREE.Vector3().subVectors(particles[i].position, particles[j].position).normalize();
                                pressureForce.add(direction.multiplyScalar(-particles[j].mass * pressure * gradient));

                                const viscosity = mu * (particles[j].mass / Math.max(particles[j].density, 0.0001));
                                const relativeVelocity = new THREE.Vector3().subVectors(particles[j].velocity, particles[i].velocity);
                                viscosityForce.add(relativeVelocity.multiplyScalar(-viscosity * gradient));
                            }
                        }
                    }
                }
            }
        }

        particles[i].forces.add(pressureForce);
        particles[i].forces.add(viscosityForce);
        particles[i].forces.add(gravityForce);
    }
}

// Function to calculate the forces acting on the particles using the naive approach
function calculateForces(particles, h, mu, g) {
    for (let i = 0; i < particles.length; i++) {
        // Initialize the forces acting on the particle
        let pressureForce = new THREE.Vector3(0, 0, 0);
        let viscosityForce = new THREE.Vector3(0, 0, 0);
        const gravityForce = new THREE.Vector3(0, -g * particles[i].mass, 0); // Gravity is constant

        // Reset forces
        particles[i].forces.set(0, 0, 0);

        // Loop over all the particles to calculate the forces acting on the current particle
        for (let j = 0; j < particles.length; j++) {
            if (i === j) continue; // Skip the particle itself

            // Calculate r as the distance between the two particles
            const r = particles[i].position.distanceTo(particles[j].position);

            if (r < h) { // Only consider neighbors within smoothing length
                // Compute the gradient of the kernel function
                const gradient = gradientKernel(r, h);

                // Compute the pressure
                const pressure = (particles[i].pressure + particles[j].pressure) / (2 * particles[j].density);
                const direction = new THREE.Vector3().subVectors(particles[i].position, particles[j].position).normalize();
                pressureForce.add(direction.multiplyScalar(-particles[j].mass * pressure * gradient));

                // Compute the viscosity
                const viscosity = mu * (particles[j].mass / particles[j].density);
                const relativeVelocity = new THREE.Vector3().subVectors(particles[j].velocity, particles[i].velocity);
                viscosityForce.add(relativeVelocity.multiplyScalar(-viscosity * gradient));
            }
        }


        particles[i].forces.add(pressureForce);
        particles[i].forces.add(viscosityForce);
        particles[i].forces.add(gravityForce);
    }
}

// Function to update the position of the particles based on the forces acting on them
function updateParticles(particles, timeStep) {
    for (let i = 0; i < particles.length; i++) {

        const acceleration = particles[i].forces.clone().divideScalar(particles[i].density);
        particles[i].velocity.add(acceleration.multiplyScalar(timeStep));
        particles[i].position.add(particles[i].velocity.clone().multiplyScalar(timeStep));

        if (!simulationParams.InstancedMesh) {
            // Update the position of the particle mesh
            particles[i].mesh.position.copy(particles[i].position);
            // Change the color of the particles based on their velocity
            const velocity = particles[i].velocity.length();
            const color = new THREE.Color(0x0000ff);
            color.setHSL((velocity + 1) / 2, 1, 0.5);
            particles[i].sphereMaterial.color = color;
        } else {
            const dummy = new THREE.Object3D();
            dummy.position.copy(particles[i].position);
            dummy.updateMatrix();
            instancedMesh.setMatrixAt(particles[i].index, dummy.matrix);
        }

        // Check for collisions with the simulation boundary
        if (particles[i].position.x < simulationParams.boxWidth / -2 || particles[i].position.x > simulationParams.boxWidth / 2) {
            particles[i].velocity.x *= -1 * simulationParams.BoxRestitution;
            particles[i].position.x = Math.max(Math.min(particles[i].position.x, simulationParams.boxWidth / 2), simulationParams.boxWidth / -2);
        }
        if (particles[i].position.y < simulationParams.boxHeight / -2 || particles[i].position.y > simulationParams.boxHeight / 2) {
            particles[i].velocity.y *= -1 * simulationParams.BoxRestitution;
            particles[i].position.y = Math.max(Math.min(particles[i].position.y, simulationParams.boxHeight / 2), simulationParams.boxHeight / -2);
        }
        if (simulationParams.ThirdDimension) {
            if (particles[i].position.z < simulationParams.boxDepth / -2 || particles[i].position.z > simulationParams.boxDepth / 2) {
                particles[i].velocity.z *= -1 * simulationParams.BoxRestitution;
                particles[i].position.z = Math.max(Math.min(particles[i].position.z, simulationParams.boxDepth / 2), simulationParams.boxDepth / -2);
            }
        }
        
    }
}

// Function to compute the reflections of the particles if enabled
function computeReflections(particles) {
    if (simulationParams.enableReflections) {
        for (let i = 0; i < particles.length; i++) {
            const particle = particles[i];
            particle.mesh.visible = false; // Hide the particle during the reflection capture
            particle.cubeCamera.position.copy(particle.position);
            particle.cubeCamera.update(renderer, scene);
            particle.mesh.visible = true; // Show the particle after the reflection capture
        }
    }
}

// Three.js rendering loop function
function animate() {

    if (simulationParams.InstancedMesh) {
        instancedMesh.instanceMatrix.needsUpdate = true;
    }

    requestAnimationFrame(animate);

    if (simulationParams.useGridOptimization) {
        calculateDensityAndPressureOptimized(
            particles,
            simulationParams.smoothingLength,
            simulationParams.referenceDensity, 
            simulationParams.gasConstant
        );
        calculateForcesOptimized(
            particles, 
            simulationParams.smoothingLength, 
            simulationParams.referenceDensity, 
            simulationParams.gravity
        );
    } else {
        calculateDensityAndPressure(
            particles,
            simulationParams.smoothingLength,
            simulationParams.referenceDensity,
            simulationParams.gasConstant
        );
        calculateForces(
            particles,
            simulationParams.smoothingLength,
            simulationParams.referenceDensity,
            simulationParams.gravity
        );
    }

    updateParticles(particles, simulationParams.timeStep);
    computeReflections(particles);

    stats.update();
    renderer.render(scene, camera);
}

// Function to rotate the camera using the mouse
function rotateCameraMouse() {
    let isMouseDown = false;
    let previousMousePosition = {
        x: 0,
        y: 0
    };

    window.addEventListener('mousedown', (e) => {
        isMouseDown = true;
        previousMousePosition = {
            x: e.clientX,
            y: e.clientY
        };
    });

    window.addEventListener('mouseup', (e) => {
        isMouseDown = false;
    });

    window.addEventListener('mousemove', (e) => {

        // if the mouse is on the dat.gui panel, don't rotate the camera
        const gui = document.querySelector('.dg.ac');
        if (gui && gui.contains(e.target)) {
            return;
        }

        if (isMouseDown) {
            const deltaMove = {
                x: e.clientX - previousMousePosition.x,
                y: e.clientY - previousMousePosition.y
            };

            camera.rotation.y += deltaMove.x * 0.01;
            camera.rotation.x += deltaMove.y * 0.01;

            previousMousePosition = {
                x: e.clientX,
                y: e.clientY
            };
        }
    });
}

// Function to move the camera using the keyboard arrows
function moveCameraArrows() {
    window.addEventListener('keydown', (e) => {
        switch (e.key) {
            case 'w':
                camera.position.z -= 0.5;
                simulationParams.zSimulation = camera.position.z;
                break;
            case 's':
                camera.position.z += 0.5;
                simulationParams.zSimulation = camera.position.z;
                break;
            case 'a':
                camera.rotation.y += 0.25;
                break;
            case 'd':
                camera.rotation.y -= 0.25;
                break;
        }
    });
}

// Function to resize the renderer when the window is resized
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}


// Setup a scene, camera, renderer
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer();

// set camera position
camera.position.z = simulationParams.zSimulation;

//set background color
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0x33393c)
document.body.appendChild(renderer.domElement);

// Add to the DOM the stats panel
const stats = new Stats();
document.body.appendChild(stats.dom);

// Add an ambient light to the scene
const ambientLight = new THREE.AmbientLight(0x404040, 0.5); // Soft white light
scene.add(ambientLight);

// Add a directional light to the scene
const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
directionalLight.position.set(1, 1, 1).normalize();
scene.add(directionalLight);

// Load the environment map
const loader = new THREE.TextureLoader();
loader.load('map.jpg', function (texture) {
    texture.mapping = THREE.EquirectangularReflectionMapping;
});

// Create a geometry for the particles and the box
let sphereGeometry, cubeGeometry;
if (simulationParams.ThirdDimension) {
    sphereGeometry = new THREE.SphereGeometry(simulationParams.particlesDimensions, 16, 16);
    cubeGeometry = new THREE.BoxGeometry(simulationParams.boxWidth, simulationParams.boxHeight, simulationParams.boxDepth);
} else {
    sphereGeometry = new THREE.CircleGeometry(simulationParams.particlesDimensions, 16);
    cubeGeometry = new THREE.BoxGeometry(simulationParams.boxWidth, simulationParams.boxHeight, 0.0);
}


// Create a grid of particles
const particles = gridOfParticles(
    simulationParams.numberOfParticles,
    simulationParams.numberOfParticles,
    simulationParams.ThirdDimension ? simulationParams.numberOfParticles : 1,
    simulationParams.boxWidth,
    simulationParams.boxHeight,
    simulationParams.boxDepth,
    simulationParams.particleSpacing
);


// Create a cube as bounding box for the simulation
let cubeMaterial = new THREE.MeshBasicMaterial({ color: 0x000000, wireframe: true });
let cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
scene.add(cube);

// Create an InstancedMesh for rendering the particles if the instancedMesh parameter is enabled
const instancedMaterial = new THREE.MeshPhysicalMaterial({ color: 0xffff00, metalness: simulationParams.particlesMetalness, roughness: simulationParams.particlesRoughness });
if (simulationParams.InstancedMesh) {
    var instancedMesh; // This is a var, so it can be accessed from the animate function (global scope)
    instancedMesh = new THREE.InstancedMesh(sphereGeometry, instancedMaterial, particles.length);
    scene.add(instancedMesh);
}

updateBoxGeometry();
window.addEventListener('resize', onWindowResize, false);
rotateCameraMouse();
moveCameraArrows();
createGui();
animate();