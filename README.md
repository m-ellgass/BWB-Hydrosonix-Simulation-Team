# BWB Hydrosonix Simulation Team

Repository created in 2026 to house and keep updated all of the relevant shear wave/simulation team code for BWB Hydrosonix. The main file for output graphs is `process_model_1d_recon.m`, which tests how accurately three different time-of-flight (TOF) methods recover wave group velocity from a simulated 1D wavefield with noise added. It loads a simulation, crops it, adds noise at a chosen SNR, band-pass filters it, then estimates group velocity three ways: NCC, envelope fit, and TTP. The simulated data is compared to a "true" shear wave speed, held in 1d_sim_cs_x.mat, which is created by `run_model_1d.m`.

## Repository Structure

  `sw1d.m`: the shear wave propagation code for 1D, outputs displacement, mesh size, and time points

  `run_model_1d.m`: the main simulation file that outputs the shear wave results, specifies simulation parameters, and calls sw1d.m to solve

      ^ this simulation file outputs '1d_sim_cs_x.mat`, where x is the simulated speed of the wave

  `process_model_1d_recon.m`: takes simulation data and processes it with several techniques (white noise, correlation, envelope fitting) and plots results. 

Main files call the following functions:

    `sim_crop.m`: to resize array for analysis

    `sd_filter_surface_wavefield.m`: applies filter to surface wave to then plot the filtered wave propagation

    `pulse_width.m`: solves for pulse characteristics and signal to noise ratio (currently commented out)

    `xcorr_td.m`: function to find correlation matrix between two curves, used in xttools_group_velocity_ncc.m

    `xttools_power_spectrum.m`: solves for power spectrum

    `xttools_group_velocity_ncc.m`: solves for group velocity with ncc kernel technique, used to calculate time of flight

    `xttools_group_velocity_env.m`: solves for group velocity with envelope fitting, used to calculate time of flight

    `xttools_group_velocity_ttp.m`: solves for group velocity with ttp technique, used to calculate time of flight
