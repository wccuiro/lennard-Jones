use rand::{Rng, seq::SliceRandom};

#[allow(dead_code)]
struct Universe {
    epsilon: f64,
    sigma: f64,
    beta: f64,
    particles: usize,
    length: f64
}

impl Universe {
    fn distance(&self, position_in: &[f64;2], position_fi: &[f64;2]) -> f64 {
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

        (x+y).sqrt()
    }

    fn potential(&self, index: usize, particle: &[f64;2], set: &Vec<[f64;2]>) -> f64 {
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

    fn pressure_particle(&self, index: usize, particle: &[f64;2], set: &Vec<[f64;2]>) -> f64 {
        let mut vir: f64 = 0.0;
        for (index_i, &i) in set.iter().enumerate() {
            let d: f64 = self.distance(&particle,&i);
            
            if (d >= self.length/2.0) | (index_i==index) {
                vir += 0.0 
            } else {
                vir += 8.0 * self.epsilon * ( 2. * (self.sigma/d).powi(12) - (self.sigma/d).powi(6))
            }
        };
        let density: f64 = (self.particles as f64)/self.length.powi(2);
        density/self.beta + vir/self.length.powi(2)
    }

    fn boundary(&self, r: &[f64;2]) -> [f64;2] {
        let mut rp: [f64;2] = [0.0, 0.0];

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

        rp
    }

    fn init_position(&self) -> Vec<[f64;2]> {
        let lim :f64 = (usize::max_value() as f64).powf(-1./48.);
        let sites :usize = ( self.length / lim ).floor() as usize + 1;
        
        let side = self.length / sites as f64;

        if sites.pow(2) <= self.particles {
            println!("Warning");
        }

        let mut initial: Vec<f64> = vec![0.0; sites * sites - self.particles];
        initial.extend(vec![1.0; self.particles]);

        let mut rng = rand::thread_rng();
        initial.shuffle(&mut rng);
        let mut r: Vec<[f64;2]> = vec![[0.0; 2]; self.particles];
        let mut k = 0;
        
        for (i , &j) in initial.iter().enumerate() {
            if j == 1.0 {
                let x = (i % sites) as f64 * side + 0.5 * side;
                let y = (i / sites) as f64 * side + 0.5 * side;
                r[k][0] = x;
                r[k][1] = y;
                k += 1;
            }
        };
        r
    }

    fn energy(&self, set: &Vec<[f64;2]>) -> Vec<f64> {
        let mut energy: Vec<f64> = vec![0.0; self.particles];
        for k in 0..self.particles {
            energy[k] = self.potential(k, &set[k], set);
        };
        energy
    }

    fn pressure(&self, set: &Vec<[f64;2]>) -> Vec<f64> {
        let mut vir: Vec<f64> = vec![0.0; self.particles];
        for k in 0..self.particles {
            vir[k] = self.pressure_particle(k, &set[k], set);
        };
        vir
    }

    fn movement(&self, set: &mut Vec<[f64;2]>, energy: &mut Vec<f64>, pressure: &mut Vec<f64>) {
        let delta: f64 = 0.25;
        
        let mut rng = rand::thread_rng();
        let i:usize = rng.gen_range(0..self.particles);

        let random_float1: f64 = rng.gen();
        let random_float2: f64 = rng.gen();
        let random_float3: f64 = rng.gen();
        
        let mut pos_final: [f64;2] = [0.0;2];
        pos_final[0] = set[i][0] + (random_float1-0.5) * delta;
        pos_final[1] = set[i][1] + (random_float2-0.5) * delta;
        
        pos_final = self.boundary(&pos_final);

        let energy_move: f64 = self.potential(i, &pos_final, set);
        let delta_energy: f64 = energy_move - energy[i];

        if random_float3 <= (-self.beta*(delta_energy)).exp() {
            set[i] = pos_final;
            energy[i] = energy_move;
            pressure[i] = self.pressure_particle(i, &pos_final, set);
        }
    }

    fn total_energy(&self, energy: &Vec<f64>) -> f64 {
        let mut total: f64 = 0.;
        for i in energy {
            total += i;
        };
        total/2.
    }

    fn total_pressure(&self, pressure: &Vec<f64>) -> f64 {
        let mut total: f64 = 0.;
        for i in pressure {
            total += i;
        }
        let density: f64 = (self.particles as f64)/self.length.powi(2);
        density/self.beta + total/(2.0 * self.length.powi(2))
    }
}
//use std::fs::File;
//use std::io::prelude::*;
//use std::path::Path;

use statrs::statistics::Statistics;

fn simulation(length: f64) -> (f64, f64, f64) {
    let temperature:f64 = 0.7; 
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
    
    let mut positions: Vec<[f64;2]> = box1.init_position();
    let mut energies: Vec<f64> = box1.energy(&positions);
    let mut pressures: Vec<f64> = box1.pressure(&positions);
    
    let mov_eq: usize = 2500000;
    let mov_dec: usize = 10;
    let mov_samples: usize = 10000;

    let mut energy: f64 = 0.0;
    let mut pressure: f64 = 0.0;

    let mut energy_samples: Vec<f64> = Vec::with_capacity(mov_samples);

    let mut pressure_samples: Vec<f64> = Vec::with_capacity(mov_samples);

    for _ in 0..mov_eq {
        box1.movement(&mut positions, &mut energies, &mut pressures);
    }

    for i in 0..(mov_dec * mov_samples) {
        box1.movement(&mut positions, &mut energies, &mut pressures);
        
        energy += box1.total_energy(&energies)/(mov_dec as f64);
        pressure += box1.total_pressure(&pressures)/(mov_dec as f64);
            
        if i % mov_dec == 0 {
            energy_samples.push(energy);
            pressure_samples.push(pressure);
                
//                match write!(file, "{}\n", energy/(mov_dec as f64)) {
//                    Err(why) => panic!("couldn't write to {}:{}", display, why),
//                    Ok(_) => (),
//                };
            //println!("{}",energy);
            energy = 0.;
            pressure = 0.;
        }
    }

    let energy_mean: f64 = (&energy_samples).mean();
    let energy_standard_deviation: f64 = (&energy_samples).std_dev();
    let pressure_mean: f64 = pressure_samples.mean();
    (energy_mean, energy_standard_deviation, pressure_mean)
}

use rayon::prelude::*;

fn main() {
    let init: f64 = 22.0 ;
    let points: usize = 14 ;
    let end: f64 = 70.0 ;

    let step: f64 = (end-init)/(points as f64);

    (0..=points).into_par_iter().for_each(|i| {
        let length: f64 = init + (i as f64)*step;
        let (energy, heat_capacity, pressure) : (f64, f64, f64) = simulation(length);
        
        let density: f64 = 500.0/length.powi(2);
        //println!("{:?}", energies)
        println!("{}, {}, {}, {}", energy, heat_capacity, pressure, density);
    });
}
