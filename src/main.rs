use rand::{Rng, seq::SliceRandom};

/// The universe class is defined here, what is the type of the variables. The value of the
/// constants can also be defined here.
#[allow(dead_code)]
struct Universe {
    epsilon: f64,
    sigma: f64,
    beta: f64,
    particles: usize,
    length: f64
}

/// Here the methods of the Universe class are implemented.
impl Universe {

    /// Computes the periodic distance betweeen two atoms in three dimensions. 
    ///
    /// $$
    /// d = min(||\vec{x} - \vec{y}||, ||(\vec{x} - \vec{y}) \pm (L,0,0)||, ||(\vec{x} - \vec{y})
    /// \pm (0,L,0)||, ||(\vec{x} - \vec{y}) \pm (0,0,L)||)
    /// $$
    ///
    /// # Example
    ///
    /// ```
    /// dis =  distance(&[2.,4.,5.], &[3.,4.,1.])
    /// ```

    fn distance(&self, position_in: &[f64;3], position_fi: &[f64;3]) -> f64 {
        let x = {
            let mut positions = vec![
            (position_in[0]-position_fi[0]).powi(2), 
            (position_in[0]-position_fi[0] - self.length ).powi(2),
            (position_in[0]-position_fi[0] + self.length ).powi(2)
            ];
            positions.sort_by(|a, b| a.partial_cmp(b).unwrap());
            positions[0]
        };
        
        let y = {
            let mut positions = vec![
            (position_in[1]-position_fi[1]).powi(2), 
            (position_in[1]-position_fi[1] - self.length ).powi(2),
            (position_in[1]-position_fi[1] + self.length ).powi(2)
            ];
            positions.sort_by(|a, b| a.partial_cmp(b).unwrap());
            positions[0]
        };

        let z = {
            let mut positions = vec![
            (position_in[2]-position_fi[2]).powi(2), 
            (position_in[2]-position_fi[2] - self.length ).powi(2),
            (position_in[2]-position_fi[2] + self.length ).powi(2)
            ];
            positions.sort_by(|a, b| a.partial_cmp(b).unwrap());
            positions[0]
        };
       

        (x+y+z).sqrt()
    }

    /// # Computes the energy of one particle in the system
    ///
    /// The energy is only computed if the periodic distance is less than L. We are using the
    /// formula. But tail corrections are not implemnted.
    ///
    /// $$
    /// E = 4 \epsilon [(\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6}]
    /// $$
    ///
    /// where $ \epsilon, \sigma $ are constants defined in the class and $r$ is the periodic
    /// distance.

    fn energy_particle(&self, index: usize, particle: &[f64;3], set: &Vec<[f64;3]>) -> f64 {
        let mut energy: f64 = 0.0;
        for (index_i, &i) in set.iter().enumerate() {
            let d: f64 = self.distance(&particle,&i);
            
            if (d >= self.length/2.0) | (index_i==index) {
                energy += 0.0 
            } else {
                energy += 4.0 * self.epsilon * ((self.sigma/d).powi(12)-(self.sigma/d).powi(6))
            }
        };
        energy
    }

    /// # Computes the virial of one particle in the system
    ///
    /// The virial is computed assuming the full curve of the Lennard Jones potential. But no tail
    /// corrections are implemented.
    ///
    /// $$
    /// Vir = 8 \epsilon [2(\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6}]
    /// $$
    /// 
    /// where $\epsilon, \sigma$ are constants defined in the class and $r$ is the periodic
    /// distance.

    fn pressure_particle(&self, index: usize, particle: &[f64;3], set: &Vec<[f64;3]>) -> f64 {
        let mut vir: f64 = 0.0;
        for (index_i, &i) in set.iter().enumerate() {
            let d: f64 = self.distance(&particle,&i);
            
            if (d >= self.length/2.0) | (index_i==index) {
                vir += 0.0 
            } else {
                vir += 8.0 * self.epsilon * ( 2. * (self.sigma/d).powi(12) - (self.sigma/d).powi(6))
            }
        };
        vir
    }

    /// # Keeping everything in the box
    ///
    /// Since in each movement is possible to move out of the box. This method applies periodic
    /// boundary conditions to the coordinates of the particles if is outside the limits.
    ///

    fn boundary(&self, r: &[f64;3]) -> [f64;3] {
        let mut rp: [f64;3] = [0.0;3];

        if r[0] > self.length {
            rp[0] = r[0] - self.length
        } else if r[0] < 0. {
            rp[0] = r[0] + self.length
        } else {
            rp[0] = r[0]
        }

        if r[1] > self.length {
            rp[1] = r[1] - self.length
        } else if r[1] < 0. {
            rp[1] = r[1] + self.length
        } else {
            rp[1] = r[1]
        }

        if r[2] > self.length {
            rp[2] = r[2] - self.length
        } else if r[2] < 0. {
            rp[2] = r[2] + self.length
        } else {
            rp[2] = r[2]
        }

       rp
    }

    /// # Initial Configuration
    ///
    /// Knowing the number of particles and the length of the box this method puts the particles in
    /// a radom position. To do that first computes the maximum ammount of particles that can be
    /// inside the box. This limit depends on the value of the highest number that is manageble by the system.
    ///
    /// $$
    /// N_{max} = (\frac{L}{k^{-1/48}})^{3}
    /// $$
    /// 
    /// where $k$ is the maximum unisigned integer the system can manage.
    ///
    /// With this number We can discretize the space creating $N_{max}$ cells and then putting a
    /// particle randomly inside them. Finally returning only the position.
    ///
    fn init_position(&self) -> Vec<[f64;3]> {
        let lim :f64 = (usize::max_value() as f64).powf(-1./48.);
        let sites :usize = ( self.length / lim ).floor() as usize + 1;
        let side = self.length / sites as f64;

        if sites.pow(3) <= self.particles {
            println!("Warning");
        }

        let mut initial: Vec<f64> = vec![0.0; sites.pow(3) - self.particles];
        initial.extend(vec![1.0; self.particles]);

        let mut rng = rand::thread_rng();
        initial.shuffle(&mut rng);
        let mut r: Vec<[f64;3]> = vec![[0.0; 3]; self.particles];
        let mut k = 0;
        
        for (i , &j) in initial.iter().enumerate() {
            if j == 1.0 {
                let z = (i / sites.pow(2)) as f64 * side + 0.5 * side;
                let layer: usize = i % sites.pow(2);
                let x = (layer % sites) as f64 * side + 0.5 * side;
                let y = (layer / sites) as f64 * side + 0.5 * side;
                r[k][0] = x;
                r[k][1] = y;
                r[k][2] = z;
                k += 1;
            }
        };
        r
    }

    /// # Total energy
    /// 
    /// Computes the total energy of the system calling the function for the energy of one
    /// particle and at the end dividing the result by 2.
    ///

    fn energy(&self, set: &Vec<[f64;3]>) -> f64 {
        let mut energy: f64 = 0.0;
        for k in 0..self.particles {
            energy += self.energy_particle(k, &set[k], set);
        };
        energy/2.0
    }

    /// # Total pressure
    /// 
    /// Computes the total pressure of the system calling the function for the virial of one
    /// particle. Then using
    ///
    /// $$
    /// P = \frac{\rho}{\beta} + \frac{1}{2 V} \sum Vir
    /// $$
    /// 
    /// where $N$ is the number of particles, $V=L^{3}$ is the volume, $\rho = N/V$ the density and $\beta$ is the inverse temperature.
    ///

   fn pressure(&self, set: &Vec<[f64;3]>) -> f64 {
        let mut vir: f64 = 0.0;
        for k in 0..self.particles {
            vir += self.pressure_particle(k, &set[k], set);
        };
        let density: f64 = (self.particles as f64)/self.length.powi(3);
        density/self.beta + vir/(2.0 * self.length.powi(3))
    }

   /// # Montecarlo Movement
   ///
   /// This method takes an initial set of positions, initial value of energy, initial value of pressure.
   /// Pick a random particle of this set of positions and propose a movement, subtract the
   /// energies and then compare it with a random number.
   ///
   /// If $x \leq e^{-\beta (\Delta E)}$ then the movement is accepted and the values of energy and
   /// pressure are updated.
   ///
    fn movement(&self, set: &mut Vec<[f64;3]>, energy: &mut f64, pressure: &mut f64) {
        let delta: f64 = 0.25;
        
        let mut rng = rand::thread_rng();
        let i:usize = rng.gen_range(0..self.particles);

        let random_float1: f64 = rng.gen();
        let random_float2: f64 = rng.gen();
        let random_float3: f64 = rng.gen();
        
        let mut pos_final: [f64;3] = [0.0;3];
        pos_final[0] = set[i][0] + (random_float1-0.5) * delta;
        pos_final[1] = set[i][1] + (random_float2-0.5) * delta;
        pos_final[2] = set[i][2] + (random_float2-0.5) * delta;
        
        pos_final = self.boundary(&pos_final);

        let energy_init: f64 = self.energy_particle(i, &set[i], set);
        let energy_move: f64 = self.energy_particle(i, &pos_final, set);
        let delta_energy: f64 = energy_move - energy_init;

        if random_float3 <= (-self.beta*(delta_energy)).exp() {
            let pressure_init: f64 = self.pressure_particle(i, &set[i], set);
            let pressure_move: f64 = self.pressure_particle(i, &pos_final, set);

            let delta_pressure: f64 = pressure_move - pressure_init;

            set[i] = pos_final;
            *energy += delta_energy;
            *pressure += delta_pressure/self.length.powi(3);
        }
    }
}
//use std::fs::File;
//use std::io::prelude::*;
//use std::path::Path;

use statrs::statistics::Statistics;

/// # Montecarlo simulation
///
/// This function run a simulation for a given length of the box. It starts initializing the struct
/// and then initializing the initial positions, energy and pressure. Then run some movements
/// called equilibrium points where no data is saved. Then takes data and saving it every mov_dec
/// amount of movements.
///
/// Finally returns the average energy, stantdard deviation of the energy (heat capacity), average
/// pressure and standard deviation of the pressure.
///
fn simulation(length: f64) -> (f64, f64, f64, f64) {
    let temperature:f64 = 1.5; 
//    let base_filename = "./iteration";

//    let filename = format!("{}{}.dat",base_filename,iteration);
//    let path = Path::new(&filename);
//    let display = path.display();

//    let mut file = match File::create(&path) {
//        Err(why) => panic!("couldn't create {}: {}", display, why),
//        Ok(file) => file,
//    };

    let box1 = Universe {
        epsilon: 1.,
        sigma: 1.,
        beta: 1./temperature,
        particles: 500,
        length,
    };
    
    let mut positions: Vec<[f64;3]> = box1.init_position();
    let mut energy: f64 = box1.energy(&positions);
    let mut pressure: f64 = box1.pressure(&positions);
    
    let mov_eq: usize = 2500000;
    let mov_dec: usize = 10;
    let mov_samples: usize = 10;

    let mut energy_ac: f64 = 0.0;
    let mut pressure_ac: f64 = 0.0;

    let mut energy_samples: Vec<f64> = Vec::with_capacity(mov_samples);

    let mut pressure_samples: Vec<f64> = Vec::with_capacity(mov_samples);

    for i in 0..mov_eq {
        box1.movement(&mut positions, &mut energy, &mut pressure);
        println!("{}, {}", i, energy)
    }

    for i in 0..(mov_dec * mov_samples) {
        box1.movement(&mut positions, &mut energy, &mut pressure);
        
        pressure_ac += pressure/(mov_dec as f64);
        energy_ac += energy/(mov_dec as f64);
            
        if i % mov_dec == 0 {
            energy_samples.push(energy_ac);
            pressure_samples.push(pressure_ac);
                
//                match write!(file, "{}\n", energy/(mov_dec as f64)) {
//                    Err(why) => panic!("couldn't write to {}:{}", display, why),
//                    Ok(_) => (),
//                };
            //println!("{}",energy);
            energy_ac = 0.;
            pressure_ac = 0.;
        }
    }

    let energy_mean: f64 = (&energy_samples).mean();
    let energy_standard_deviation: f64 = (&energy_samples).std_dev();
    let pressure_mean: f64 = (&pressure_samples).mean();
    let pressure_std: f64 = (&pressure_samples).std_dev();
    (energy_mean, energy_standard_deviation, pressure_mean, pressure_std)
}

use rayon::prelude::*;

/// # Main function
///
/// It defines the initial and final density where the simulation will run. Also the amount of
/// points that will be taken in this range. Finally prints 
/// (density, energy, energy standard deviation, pressure, pressure standard deviation)
///
/// It run every density in parallel.

fn main() {
    let init_density: f64 = 0.3 ;
    let points: usize = 1 ; // In reality the total number of points is points +1
    let end_density: f64 = 0.3 ;

    let step: f64 = (end_density-init_density)/(points as f64);

    (0..points).into_par_iter().for_each(|i| {
        let density: f64 =  init_density + (i as f64) * step;
        let length: f64 = (500.0/density).powf(1./3.);
        let (energy, heat_capacity, pressure, std_pressure) : (f64, f64, f64, f64) = simulation(length);
        
//        println!("{:?}", energies)
//        println!("{}, {}, {}, {}, {}", density, energy, heat_capacity, pressure, std_pressure);
//        println!("{}", length)
    });
}
