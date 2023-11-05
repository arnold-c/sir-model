# %%
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from plotnine import *


# %%
def seirmod(t, y, beta, mu, sigma, gamma, N):
    # Unpack states
    S, E, I, R = y

    # Define equations
    dS = mu * (N - S) - beta * S * I / N
    dE = beta * S * I / N - sigma * E
    dI = sigma * E - (mu + gamma) * I
    dR = gamma * I - mu * R

    # Return list of gradients
    return dS, dE, dI, dR


# %%
tmin = 0
tmax = 26
tstep = 1 / 10
times = np.arange(tmin, tmax, tstep)

beta = 2
mu = 0
sigma = 1
gamma = 1 / 2
N = 1
parms = (beta, mu, sigma, gamma, N)

S0 = 0.999
E0 = 0
I0 = 0.001
R0 = 0
start = (S0, E0, I0, R0)

# %%
out = solve_ivp(seirmod, [tmin, tmax], np.array(start), args=parms, t_eval=times)

# %%
out_df = (
    pd.DataFrame(out.y).transpose().rename(columns={0: "S", 1: "E", 2: "I", 3: "R"})
)
out_df["time"] = out.t
out_df = out_df.melt(id_vars="time", value_vars=["S", "E", "I", "R"]).rename(
    columns={"variable": "state", "value": "number"}
)

# %%
theme_set(theme_minimal())

(
    ggplot(out_df, aes(x="time", y="number", color="state"))
    + geom_line(size=2)
    + labs(x="Time", y="Number", color="State")
)

# %%
# Candidate values for R0 and beta
R0 = np.linspace(0.1, 5, 50)
betas = R0 * 1 / 2

# %%
solve_ivp(seirmod, [tmin, 1e5], start, args=parms).y[2, -1]

# %%
final_size_df = pd.DataFrame({"R0": R0, "final_size": np.zeros(len(R0))})

for index, beta in enumerate(betas):
    p = (beta, mu, sigma, gamma, N)
    final_size_df.final_size[index] = solve_ivp(seirmod, [tmin, 1e5], start, args=p).y[
        2, -1
    ]

# %%
(
    ggplot(final_size_df, aes(x="R0", y="final_size"))
    + geom_line(size=2)
    + labs(x="R0", y="Final size")
)
