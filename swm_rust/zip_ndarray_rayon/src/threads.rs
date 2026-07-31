use std::collections::BTreeSet;
use std::io::Write;

// supplied by Rory Kelly

/// How to bind Rayon's worker threads to CPUs.
#[derive(Clone, Copy)]
enum Affinity {
    /// Leave threads unbound; the OS scheduler moves them freely (default).
    None,
    /// Pin one thread to each logical CPU (includes SMT/hyperthread siblings).
    Logical,
    /// Pin one thread per physical core (skips SMT/hyperthread siblings).
    Physical,
}

/// Parsed command-line options.
struct Config {
    affinity: Affinity,
    /// Explicit pool size from `--threads`; `None` means use the default.
    threads: Option<usize>,
}

const HELP: &str = "\
rust-threads — a small Rayon threading demo

Usage: rust-threads [OPTIONS]

Options:
  --affinity=<none|logical|physical>
        Bind worker threads to CPUs (default: none)
          none      leave threads unbound, let the OS schedule them
          logical   pin one thread per logical CPU (incl. SMT siblings)
          physical  pin one thread per physical core (skip SMT siblings)
  --threads=<N>
        Number of threads in the pool (default: one per available CPU core)
  -h, --help
        Show this help and exit";

pub fn threads() {
    let config = parse_args();

    let cpus = std::thread::available_parallelism().unwrap();
    // "CPUs" sorts before "Hello", so a plain `sort` keeps this line on top.
    println!("CPUs visible: {cpus}");

    // The CPUs to pin to (empty for `none`), and the default pool size for the
    // chosen affinity mode (one thread per CPU/core that mode targets).
    let pinned = cores_to_pin(config.affinity);
    let default_threads = if pinned.is_empty() {
        cpus.get()
    } else {
        pinned.len()
    };
    let num_threads = config.threads.unwrap_or(default_threads);

    let mut builder = rayon::ThreadPoolBuilder::new().num_threads(num_threads);
    if !pinned.is_empty() {
        // Rayon calls this once on each worker as it starts up; we use it to
        // bind that worker to a CPU. If there are more threads than CPUs we
        // wrap around so several threads share a core.
        builder = builder.start_handler(move |index| {
            core_affinity::set_for_current(pinned[index % pinned.len()]);
        });
    }
    builder
        .build_global()
        .expect("failed to configure Rayon thread pool");

    let mode = match config.affinity {
        Affinity::None => "none (unbound, OS-scheduled)",
        Affinity::Logical => "logical (one thread per logical CPU)",
        Affinity::Physical => "physical (one thread per physical core)",
    };
    let num_threads = rayon::current_num_threads();
    println!("Affinity mode: {mode}");
    println!("Running on {num_threads} rayon threads\n");

    // Flush the header to the OS before the parallel region starts, so it can
    // never interleave with the per-thread output below (e.g. when piped).
    std::io::stdout().flush().unwrap();

    // `broadcast` runs the closure exactly once on every thread in the pool,
    // so each worker reports itself once — no work-stealing surprises.
    rayon::broadcast(|ctx| {
        let id = ctx.index();
        // The CPU this thread is running on right now. With pinning it stays
        // fixed; unbound, it's just wherever the OS placed it this instant.
        // let cpu = unsafe { libc::sched_getcpu() };
        let cpu = -1;
        println!("Hello from rayon: thread {id}, CPU id {cpu}");
    });
}

/// Parse the command line. Exits the process on `--help` or a bad argument.
fn parse_args() -> Config {
    let mut config = Config {
        affinity: Affinity::None,
        threads: None,
    };

    for arg in std::env::args().skip(1) {
        if arg == "-h" || arg == "--help" {
            println!("{HELP}");
            std::process::exit(0);
        } else if let Some(value) = arg.strip_prefix("--affinity=") {
            config.affinity = match value {
                "none" => Affinity::None,
                "logical" => Affinity::Logical,
                "physical" => Affinity::Physical,
                other => fail(&format!("unknown --affinity value '{other}'")),
            };
        } else if let Some(value) = arg.strip_prefix("--threads=") {
            config.threads = match value.parse::<usize>() {
                Ok(0) | Err(_) => fail(&format!("--threads must be a positive integer, got '{value}'")),
                Ok(n) => Some(n),
            };
        } else {
            fail(&format!("unexpected argument '{arg}'"));
        }
    }
    config
}

/// Print an error plus a usage hint and exit with a non-zero status.
fn fail(message: &str) -> ! {
    eprintln!("error: {message}");
    eprintln!("try '--help' for usage");
    std::process::exit(1);
}

/// Decide which logical CPUs the worker threads should be pinned to.
fn cores_to_pin(affinity: Affinity) -> Vec<core_affinity::CoreId> {
    // `get_core_ids` returns every logical CPU available to this process,
    // already respecting any cpuset/affinity limits the OS imposes on us.
    let all = core_affinity::get_core_ids().unwrap_or_default();
    match affinity {
        Affinity::None => Vec::new(),
        Affinity::Logical => all,
        Affinity::Physical => {
            // Keep the first logical CPU we see for each physical core, so SMT
            // siblings collapse down to a single thread per core.
            let mut seen = BTreeSet::new();
            all.into_iter()
                .filter(|core| seen.insert(physical_core_key(core.id)))
                .collect()
        }
    }
}

/// Identify the physical core a logical CPU belongs to. SMT siblings share the
/// same key. On Linux this comes from sysfs; if that's unavailable we fall back
/// to treating every logical CPU as its own core (so physical == logical).
fn physical_core_key(cpu: usize) -> (i64, i64) {
    let read = |field: &str| -> i64 {
        std::fs::read_to_string(format!("/sys/devices/system/cpu/cpu{cpu}/topology/{field}"))
            .ok()
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(cpu as i64)
    };
    (read("physical_package_id"), read("core_id"))
}

