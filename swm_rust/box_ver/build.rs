fn main() {
    // Enforce exactly one feature is active
    let box_on     = cfg!(feature = "box");
    let vec_on     = cfg!(feature = "vec");
    let ndarray_on = cfg!(feature = "ndarray");
    let mdarray_on = cfg!(feature = "mdarray");

    let count = [box_on, vec_on, ndarray_on, mdarray_on].iter().filter(|&&x| x).count();
    if count > 1 {
        panic!("Only one container feature may be enabled at a time.");
    }
    if count == 0 {
        println!("cargo:warning=No feature specified. Defaulting to box. Use --features flag to select a different type.");
        println!("cargo:rustc-cfg=feature=\"box\"");
    }
}
