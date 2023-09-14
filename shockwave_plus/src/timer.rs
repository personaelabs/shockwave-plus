use ark_std::{end_timer, perf_trace::inner::TimerInfo as ArkTimerInfo, start_timer};
#[cfg(target_arch = "wasm32")]
use web_sys;

pub struct WebTimerInfo {
    label: &'static str,
}

pub enum TimerInfo {
    Web(WebTimerInfo),
    Cmd(ArkTimerInfo),
}

pub fn timer_start(label: &'static str) -> TimerInfo {
    #[cfg(target_arch = "wasm32")]
    {
        web_sys::console::time_with_label(label);
        TimerInfo::Cmd(WebTimerInfo { label })
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        TimerInfo::Cmd(start_timer!(|| label))
    }
}

pub fn timer_end(timer: TimerInfo) {
    match timer {
        TimerInfo::Web(t) => {
            web_sys::console::time_end_with_label(t.label);
        }
        TimerInfo::Cmd(t) => {
            end_timer!(t);
        }
    };
}
