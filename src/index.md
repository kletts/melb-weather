---
toc: false
---

# Melbourne Weather

Analysis of the daily seaosonal cycles of weather in Melbourne, using a multivariate Kalman Filter model. This tool illustrates the use of the R library `dlm` for state space modelling and the `d3` library for data visualisation. 

## Data loader

Read the data file from the data loader:

```js echo
const data = await FileAttachment("data/bom.csv").csv({ typed: true })
```

The 'data frame' is imported as a row based array: 

```js echo 
display(data)
```

## Methodology 

The data loader `data/bom.csv.R` implements a multivariate Kalman Filter signal over 4 weather stations 
around Melbourne. [see](https://darrenjw.github.io/time-series/ssm.html)


```tex
\begin{aligned}
\text{Observations:} \quad x_{i,t} &= F \cdot \theta_{t} + v_t \qquad &v_{t} \sim N(0, V) \\
\text{Transitions:} \quad \theta_{t}  &=  G \cdot \theta_{t-1} + w_t \qquad &w_t \sim N(0, W)
\end{aligned}
```

As a multivariate observation model, the measurement matrix F is given as the sum of a local linear trend and a latent cyclical component, the later specified as a trignometric series.

```tex 
F = \begin{bmatrix}
1 & 1 & 0 \\
1 & 1 & 0 \\
1 & 1 & 0 \\
1 & 1 & 0 \\
\end{bmatrix}
```

The updating of the 3 latent series is given by the matrix G, the trigonometric representation of the cyclical component combines a sine and cosine components of period 48 of period 48 (each half hour):

```tex
G = \begin{bmatrix}
1 & 0 & 0\\
0 & \cos \omega_j  & \sin \omega_j  \\
0 & -\sin \omega_j  & \cos \omega_j  \\
\end{bmatrix}
```

where ${tex`\omega_j = 2\pi / 48`}. 

<details>
<summary>R implementation code</summary>

```r
# The model representation as matrices
fn <- function(parm) { 
  # Optimal parameter update function
  mod <- dlm(
    FF=matrix(c(1,1,1,1,1,1,1,1,0,0,0,0), ncol=3, byrow=FALSE), 
    V=diag(exp(rep(parm[1], 4))), 
    GG=(dlmModPoly(1) + dlmModTrig(s=24*2, q=1))$GG, 
    W=diag(exp(parm[2:4])), 
    m0=c(0,0,0), 
    C0=diag(c(1,1,1)))
  return(mod) } 
# Finding the parameters of the unknown coefficients 
fit <- dlmMLE(data.m, parm=log(c(1,1,1,1)), build = fn, hessian = TRUE)
# Fitting the model to the data to measure the latent states
fitf <- dlmFilter(data.m, fn(fit$par)) 
```

</details>
<br>

## Common trend signal

Compares the observed weather by station with the common trend (removing the cycling component). The cyclical component stronest for temperature, is weakly statistically significant for wind speed, but there is no systematic cyclical component in the humidity and pressure. 

```js echo
const datas = view(Inputs.select(d3.group(data, d => d.metric)));
const inclci = view(Inputs.toggle({ label: "Include 95% CI" }));
```

```js echo
resize(width => { return Plot.plot({ 
  width, 
  marginBottom: 40,
  title: `3 days to ${d3.timeFormat("%-d %B %H:%M")(data[data.length-1].obstime)}`,
  color: { legend: true, 
    domain: ["Avalon", "Essendon", "Frankston", "Viewbank", "Trend"], 
    range: ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "black"] }, 
  y: { label: datas[0].metric }, 
  x: { type: `time`, label: `Local Time` },
  inset: 10,
  marks: [ 
    Plot.frame(), 
    inclci ? Plot.areaY(datas, { 
      x: "obstime", 
      y1: d => d.Trend + 1.96*d.TrendVar, 
      y2: d => d.Trend - 1.96*d.TrendVar, 
      fill: "lightgray", 
      opacity: 0.6 }) : null, 
    ["Avalon", "Essendon", "Frankston", "Viewbank"].map(
      d => Plot.lineY(datas, { 
        x: "obstime", 
        y: d, 
        stroke: () => d 
        })
      ), 
    Plot.lineY(datas, { 
      x: "obstime", 
      y: "Trend",  stroke: () => "Trend" })
  ]
  }); 
})
```

The latent cyclical component of the weather over the last 3 days is shown below, it shows the inherent daily variation from day to night in the different weather metrics, averaged across the 4 observation stations. 

```js 
resize(width => { return Plot.plot({ 
  width, 
  marginBottom: 40,
  title: `3 days to ${d3.timeFormat("%-d %B %H:%M")(data[data.length-1].obstime)}`,
  color: { legend: true }, 
  y: { label: datas[0].metric }, 
  x: { type: `time`, label: `Local Time` },
  inset: 10,
  marks: [ 
    Plot.frame(), 
    inclci ? Plot.areaY(datas, { 
      x: "obstime", 
      y1: d => d.Cycle + 1.96*Math.sqrt(d.CycleVar), 
      y2: d => d.Cycle - 1.96*Math.sqrt(d.CycleVar), 
      fill: "lightgray", 
      opacity: 0.6 }) : null, 
    Plot.lineY(datas, { 
      x: "obstime", 
      y: "Cycle" })
  ]
  }); 
})
```
