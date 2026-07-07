use cuda_std::prelude::*;

// Input/output type shared with the `rustc-cuda-basic` crate.
pub type T = f32;

#[kernel]
#[allow(improper_ctypes_definitions)]
pub unsafe fn add(a: &[T], b: &[T], c: *mut T) {
    let i = thread::index_1d() as usize;
    if i < a.len() {
        let elem = unsafe { &mut *c.add(i) };
        *elem = a[i] + b[i];
    }
}
