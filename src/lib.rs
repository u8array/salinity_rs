#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate alloc;

pub mod adapters;
pub mod chemistry;
#[cfg(feature = "std")]
pub mod error;
pub mod models;
pub mod salinity;

pub use crate::adapters::teos10::sa_from_sp;
#[cfg(feature = "std")]
pub use crate::error::AppError;
pub use crate::models::{Assumptions, Inputs};
pub use crate::salinity::calculator::{
    CalcResult, Components, DetailedResult, calc_salinity_sp_iterative, calc_salinity_sp_teos10,
    compute_summary, rho_from_sp, specific_gravity,
};
