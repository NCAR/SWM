
// Below is a dummy variable to test stuff. 
pub fn rand_util(in_val: i32) -> f64 {
    let mut x = 0.0;
    for i in 0..in_val {
        x += (i as f64).sin() * (i as f64).cos();
    }
    // Default we return the last line that has no semicolon, but we can also use the return keyword if needed.
    x
}