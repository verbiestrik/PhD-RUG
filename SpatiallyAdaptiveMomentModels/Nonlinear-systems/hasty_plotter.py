import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


OUTPUT_DIR = "Data-processing/Output"
FIG_DIR = "Data-processing/Figures"

SUMMARY_FILE = os.path.join(OUTPUT_DIR, "recharge_summary.csv")
HISTORY_FILE = os.path.join(OUTPUT_DIR, "recharge_history.csv")

os.makedirs(FIG_DIR, exist_ok=True)


def nearest_times(df, target_times):
    available = np.array(sorted(df["time"].unique()))
    selected = []
    for t in target_times:
        idx = np.argmin(np.abs(available - t))
        selected.append(available[idx])
    return sorted(set(selected))


def plot_mean_height(summary_df):
    plt.figure(figsize=(6, 4))
    plt.plot(summary_df["time"], summary_df["mean_h"])
    plt.xlabel("Time")
    plt.ylabel("Mean height")
    plt.title("Mean water height over time")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "mean_height_vs_time.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "mean_height_vs_time.pdf"))
    plt.close()


def plot_mean_velocity(summary_df):
    plt.figure(figsize=(6, 4))
    plt.plot(summary_df["time"], summary_df["mean_u_m"])
    plt.xlabel("Time")
    plt.ylabel("Mean velocity")
    plt.title("Mean velocity over time")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "mean_velocity_vs_time.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "mean_velocity_vs_time.pdf"))
    plt.close()


def plot_mean_a1(summary_df):
    plt.figure(figsize=(6, 4))
    plt.plot(summary_df["time"], summary_df["mean_a1"])
    plt.xlabel("Time")
    plt.ylabel(r"Mean $\alpha_1$")
    plt.title(r"Mean $\alpha_1$ over time")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "mean_a1_vs_time.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "mean_a1_vs_time.pdf"))
    plt.close()


def plot_height_snapshots(history_df, target_times):
    chosen_times = nearest_times(history_df, target_times)

    plt.figure(figsize=(6, 4))
    for t in chosen_times:
        snap = history_df[np.isclose(history_df["time"], t)]
        plt.plot(snap["x"], snap["h"], label=f"t = {t:.3f}")

    plt.xlabel("x")
    plt.ylabel("h")
    plt.title("Height profiles at selected times")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "height_snapshots.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "height_snapshots.pdf"))
    plt.close()


def plot_velocity_snapshots(history_df, target_times):
    chosen_times = nearest_times(history_df, target_times)

    plt.figure(figsize=(6, 4))
    for t in chosen_times:
        snap = history_df[np.isclose(history_df["time"], t)]
        plt.plot(snap["x"], snap["u_m"], label=f"t = {t:.3f}")

    plt.xlabel("x")
    plt.ylabel(r"$u_m$")
    plt.title("Velocity profiles at selected times")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "velocity_snapshots.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "velocity_snapshots.pdf"))
    plt.close()


def plot_a1_snapshots(history_df, target_times):
    chosen_times = nearest_times(history_df, target_times)

    plt.figure(figsize=(6, 4))
    for t in chosen_times:
        snap = history_df[np.isclose(history_df["time"], t)]
        plt.plot(snap["x"], snap["a1"], label=f"t = {t:.3f}")

    plt.xlabel("x")
    plt.ylabel(r"$\alpha_1$")
    plt.title(r"$\alpha_1$ profiles at selected times")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "a1_snapshots.png"), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, "a1_snapshots.pdf"))
    plt.close()


def main():
    summary_df = pd.read_csv(SUMMARY_FILE)
    history_df = pd.read_csv(HISTORY_FILE)

    print("Summary columns:", list(summary_df.columns))
    print("History columns:", list(history_df.columns))

    plot_mean_height(summary_df)
    plot_mean_velocity(summary_df)
    plot_mean_a1(summary_df)

    t_end = history_df["time"].max()
    target_times = [0.0, 0.25 * t_end, 0.5 * t_end, 0.75 * t_end, t_end]

    plot_height_snapshots(history_df, target_times)
    plot_velocity_snapshots(history_df, target_times)
    plot_a1_snapshots(history_df, target_times)

    print(f"Plots saved in: {FIG_DIR}")


if __name__ == "__main__":
    main()