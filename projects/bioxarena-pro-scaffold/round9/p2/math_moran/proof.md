# Exact consensus time of the neutral Moran process on the star graph

This note proves Conjecture (9) of Brewster, Svoboda, Roscow, Chatterjee, Tkadlec, and Nowak (*PNAS Nexus*, 2025; `arXiv:2503.09841v2`), establishing that under neutral birth-death (`bd`) updating on the star graph $S_n$ with $N = n + 1$ vertices ($n \ge 1$ leaves, no self-loops), starting from $N$ pairwise distinct types, the expected number of discrete steps until homogeneity is exactly

$$T_N^{\mathrm{bd}} = n^3 - n^2 + n H_n,$$

where $H_n = \sum_{j=1}^n \frac{1}{j}$ is the $n$-th harmonic number. The same coalescing-lineage argument solves arbitrary leaf-symmetric neutral replacement rules on $S_n$, yielding the exact death-birth (`db`) consensus time

$$T_N^{\mathrm{db}} = n H_n + n - 1$$

and, more generally, the exact $\lambda$-mixed consensus time for every $\lambda \in [0, 1]$.

---

## 1. Exact model specification

Following Brewster et al. (2025, Methods and SI Section 1):

* **Graph.** Let $n \ge 1$ be an integer and $N = n + 1 \ge 2$. Let $S_n = (V, E)$ be the simple undirected star graph on vertex set $V = \{0, 1, \dots, n\}$ with edge set $E = \{\{0, i\} : i \in [n]\}$, where $[n] := \{1, \dots, n\}$. Vertex $0$ is the center ($\deg(0) = n$) and vertices $1, \dots, n$ are the leaves ($\deg(i) = 1$). There are no self-loops and no multiple edges.
* **Initial state.** At time $t = 0$, every vertex holds a distinct type: $X_0(u) \neq X_0(v)$ for all distinct $u, v \in V$.
* **Time unit and update rules.** Each discrete step $t \to t+1$ selects one ordered directed edge $(B_{t+1} \to D_{t+1})$ with $\{B_{t+1}, D_{t+1}\} \in E$, sets $X_{t+1}(D_{t+1}) = X_t(B_{t+1})$, and leaves $X_{t+1}(u) = X_t(u)$ for all $u \neq D_{t+1}$. Under neutral fitness:
  * **birth-death (`bd`):** pick $B_{t+1} \in V$ uniformly ($1/N$), then pick $D_{t+1} \in \Gamma(B_{t+1})$ uniformly ($1/\deg(B_{t+1})$). Thus
    $$p_{\mathrm{in}}^{\mathrm{bd}} := \mathbb{P}(B=0, D=i) = \frac{1}{Nn}, \qquad p_{\mathrm{out}}^{\mathrm{bd}} := \mathbb{P}(B=i, D=0) = \frac{1}{N} \qquad (i \in [n]).$$
  * **death-birth (`db`):** pick $D_{t+1} \in V$ uniformly ($1/N$), then pick $B_{t+1} \in \Gamma(D_{t+1})$ uniformly ($1/\deg(D_{t+1})$). Thus
    $$p_{\mathrm{in}}^{\mathrm{db}} := \mathbb{P}(B=0, D=i) = \frac{1}{N}, \qquad p_{\mathrm{out}}^{\mathrm{db}} := \mathbb{P}(B=i, D=0) = \frac{1}{Nn} \qquad (i \in [n]).$$
  * **$\lambda$-mixed ($\lambda \in [0, 1]$):** with probability $\lambda$ perform a `bd` step and with probability $1-\lambda$ perform a `db` step, giving
    $$p_{\mathrm{in}}^{(\lambda)} = \frac{\lambda + n(1-\lambda)}{Nn}, \qquad p_{\mathrm{out}}^{(\lambda)} = \frac{n\lambda + (1-\lambda)}{Nn}.$$

More generally, let $(p_{\mathrm{in}}, p_{\mathrm{out}})$ be any pair with $p_{\mathrm{in}} > 0$, $p_{\mathrm{out}} \ge 0$, and $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$, defining an i.i.d. sequence of edge-replacement events with $\mathbb{P}(0 \to i) = p_{\mathrm{in}}$ and $\mathbb{P}(i \to 0) = p_{\mathrm{out}}$ for every $i \in [n]$.

* **Consensus (absorption) time.** Define
  $$\tau_{\mathrm{cons}} := \inf\{t \ge 0 : |X_t(V)| = 1\}, \qquad T_N := \mathbb{E}[\tau_{\mathrm{cons}}].$$

---

## 2. Neutral duality to a coalescing random walk

Let $\phi_{(B, D)} : V \to V$ be the random parent-pointer map associated with replacement $(B \to D)$:
$$\phi_{(B, D)}(u) := \begin{cases} B & \text{if } u = D, \\ u & \text{if } u \neq D. \end{cases}$$
Let $\phi_t := \phi_{(B_{t+1}, D_{t+1})}$ for $t \ge 0$. Then $(\phi_t)_{t \ge 0}$ are i.i.d. random maps $V \to V$, and by definition of the update,
$$X_{t+1}(u) = X_t(\phi_t(u)) \qquad (u \in V,\ t \ge 0).$$
By induction on $t \ge 0$,
$$X_t(u) = X_0\!\left((\phi_0 \circ \phi_1 \circ \cdots \circ \phi_{t-1})(u)\right) \qquad (u \in V),$$
where the empty composition at $t=0$ is the identity $\mathrm{id}_V$.

**Lemma 1 (Exact distributional duality).** *Define the backward coalescing random walk $(S_m)_{m \ge 0}$ on nonempty subsets of $V$ by $S_0 = V$ and $S_{m+1} = \phi_m(S_m)$ for $m \ge 0$, and let $\tau_{\mathrm{coal}} := \inf\{m \ge 0 : |S_m| = 1\}$. Then $\tau_{\mathrm{cons}} \stackrel{d}{=} \tau_{\mathrm{coal}}$, and hence $T_N = \mathbb{E}[\tau_{\mathrm{coal}}]$.*

*Proof.* Because $X_0$ is injective on $V$, we have $|X_t(V)| = 1$ if and only if $|(\phi_0 \circ \phi_1 \circ \cdots \circ \phi_{t-1})(V)| = 1$. Moreover, both $\{|X_t(V)| = 1\}$ and $\{|S_t| = 1\}$ are absorbing events (if $|A| = 1$, then $|\phi_t(A)| = 1$). Since $\phi_0, \dots, \phi_{t-1}$ are i.i.d., reversing their order preserves their joint distribution:
$$(\phi_0, \phi_1, \dots, \phi_{t-1}) \stackrel{d}{=} (\phi_{t-1}, \dots, \phi_1, \phi_0).$$
Consequently, for every integer $t \ge 0$,
$$\mathbb{P}(\tau_{\mathrm{cons}} \le t) = \mathbb{P}\!\left(|(\phi_0 \circ \cdots \circ \phi_{t-1})(V)| = 1\right) = \mathbb{P}\!\left(|(\phi_{t-1} \circ \cdots \circ \phi_0)(V)| = 1\right) = \mathbb{P}(|S_t| = 1) = \mathbb{P}(\tau_{\mathrm{coal}} \le t).$$
Therefore $\tau_{\mathrm{cons}}$ and $\tau_{\mathrm{coal}}$ have the same distribution on $\mathbb{Z}_{\ge 0} \cup \{\infty\}$, so $\mathbb{E}[\tau_{\mathrm{cons}}] = \mathbb{E}[\tau_{\mathrm{coal}}]$. $\square$

---

## 3. Exact state compression on $S_n$

For any nonempty subset $S \subseteq V$, write
$$\pi(S) := (k, c), \qquad k := |S \cap [n]| \in \{0, \dots, n\}, \qquad c := \mathbf{1}[0 \in S] \in \{0, 1\}.$$
The state space of $(k, c)$ is $\mathcal{K} := (\{0, \dots, n\} \times \{0, 1\}) \setminus \{(0, 0)\}$. Notice $|S| = k + c$, so the absorbing set $\{|S| = 1\}$ corresponds exactly to the two states $(0, 1)$ and $(1, 0)$, while $S_0 = V$ maps to $(n, 1)$. The remaining $2n - 1$ states, namely $\{(k, 1) : 1 \le k \le n\} \cup \{(k, 0) : 2 \le k \le n\}$, are non-absorbing.

**Lemma 2 (Strong lumpability and finite absorption).** *Under $\phi_m$ with law $(p_{\mathrm{in}}, p_{\mathrm{out}})$ satisfying $p_{\mathrm{in}} > 0$, $p_{\mathrm{out}} \ge 0$, and $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$, the projected process $Y_m := \pi(S_m)$ is a time-homogeneous Markov chain on $\mathcal{K}$ with absorbing states $\{(0, 1), (1, 0)\}$ and the following transition probabilities from every non-absorbing state:*

1. *From $(k, 1)$ with $1 \le k \le n$:*
   * *to $(k-1, 1)$ with probability $k p_{\mathrm{in}}$ (edge $0 \to i$ with $i \in S \cap [n]$; $i$ coalesces with $0$);*
   * *to $(k, 0)$ with probability $k p_{\mathrm{out}}$ (edge $i \to 0$ with $i \in S \cap [n]$; $0$ coalesces with $i$);*
   * *to $(k+1, 0)$ with probability $(n-k) p_{\mathrm{out}}$ (edge $i \to 0$ with $i \in [n] \setminus S$; $0$ moves to $i$);*
   * *stay at $(k, 1)$ with probability $(n-k) p_{\mathrm{in}}$ (edge $0 \to i$ with $i \in [n] \setminus S$; no effect).*

2. *From $(k, 0)$ with $2 \le k \le n$:*
   * *to $(k-1, 1)$ with probability $k p_{\mathrm{in}}$ (edge $0 \to i$ with $i \in S \cap [n]$; $i$ moves to empty $0$);*
   * *stay at $(k, 0)$ with probability $1 - k p_{\mathrm{in}}$ (all other edges).*

*Moreover, absorption occurs almost surely with finite expectation from every state, so the expected hitting times $t(k, c) := \mathbb{E}_{(k, c)}[\inf\{m \ge 0 : |S_m| = 1\}]$ are the unique finite solution to the first-step linear system with boundary conditions $t(0, 1) = t(1, 0) = 0$.*

*Proof.* Take any $S \subseteq V$ with $\pi(S) = (k, c)$, and put $L := S \cap [n]$ ($|L| = k$). Enumerating all $2n$ directed edges $(B \to D)$:
* If $c = 1$ ($0 \in S$): firing $(0 \to i)$ gives $\phi(S) = S \setminus \{i\}$ for $i \in L$ (so $\pi(\phi(S)) = (k-1, 1)$, total probability $k p_{\mathrm{in}}$) and $\phi(S) = S$ for $i \in [n] \setminus L$ (probability $(n-k)p_{\mathrm{in}}$). Firing $(i \to 0)$ moves the lineage at $0$ to $i$, giving $\phi(S) = S \setminus \{0\}$ for $i \in L$ (state $(k, 0)$, probability $k p_{\mathrm{out}}$) and $\phi(S) = (S \setminus \{0\}) \cup \{i\}$ for $i \in [n] \setminus L$ (state $(k+1, 0)$, probability $(n-k)p_{\mathrm{out}}$). Note that when $k=n$, the $(k+1, 0)$ probability is $0$. These probabilities sum to $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$ and depend only on $(k, 1)$.
* If $c = 0$ ($0 \notin S$, $k \ge 2$): firing $(0 \to i)$ gives $\phi(S) = (S \setminus \{i\}) \cup \{0\}$ for $i \in L$ (state $(k-1, 1)$, probability $k p_{\mathrm{in}}$) and $\phi(S) = S$ for $i \in [n] \setminus L$. Firing $(i \to 0)$ leaves $S$ unchanged because $0 \notin S$. Again these probabilities sum to $1$ and depend only on $(k, 0)$.

For finite absorption, observe that $p_{\mathrm{in}} > 0$ implies that from $(k, 0)$ ($k \ge 2$) the chain transitions to $(k-1, 1)$ with probability $k p_{\mathrm{in}} \ge 2 p_{\mathrm{in}} > 0$, and from $(k, 1)$ ($k \ge 1$) the chain transitions to $(k-1, 1)$ with probability $k p_{\mathrm{in}} \ge p_{\mathrm{in}} > 0$. Thus from every non-absorbing state there is a path of length at most $n$ to $(0, 1)$ with probability at least $p_{\mathrm{in}}^n > 0$. Hence $\tau_{\mathrm{coal}}$ is stochastically dominated by $n$ times a geometric random variable of parameter $p_{\mathrm{in}}^n$, so $\mathbb{E}_{(k,c)}[\tau_{\mathrm{coal}}] < \infty$ for all $(k, c) \in \mathcal{K}$, and the harmonic/first-step linear system $(I - P_{\mathrm{trans}}) t = \mathbf{1}$ on the $2n - 1$ transient states has a unique solution (since $P_{\mathrm{trans}}^n$ has spectral radius strictly less than $1$). $\square$

---

## 4. Exact solution of the recurrence

**Theorem 1 (Master star consensus-time formula).** *Let $n \ge 1$, $N = n + 1$, $p_{\mathrm{in}} > 0$, $p_{\mathrm{out}} \ge 0$ with $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$, and set $r := p_{\mathrm{out}} / p_{\mathrm{in}}$. Then for every $k \in \{1, \dots, n\}$,*
$$t(k, 1) = n H_k + \frac{n r (n k - 1)}{k + 1},$$
*and for every $k \in \{2, \dots, n\}$,*
$$t(k, 0) = n H_{k-1} + \frac{n + n^2 r (k - 1)}{k}.$$
*In particular, starting from all $N$ vertices having pairwise distinct types,*
$$T_N(p_{\mathrm{in}}, p_{\mathrm{out}}) = t(n, 1) = n H_n + n(n - 1) \frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}.$$

*Proof.* By Lemma 2, $t(k, c)$ is uniquely determined by the boundary conditions $t(0, 1) = 0$, $t(1, 0) = 0$ together with the linear equations:
1. For $2 \le k \le n$:
   $$t(k, 0) = 1 + \bigl(1 - k p_{\mathrm{in}}\bigr) t(k, 0) + k p_{\mathrm{in}}\, t(k - 1, 1),$$
   which simplifies immediately to
   $$t(k, 0) = \frac{1}{k p_{\mathrm{in}}} + t(k - 1, 1) \qquad (2 \le k \le n). \tag{1}$$
2. For $k = 1$ (note $(1, 0)$ and $(0, 1)$ are absorbing, so $t(0, 1) = t(1, 0) = 0$):
   $$\bigl(1 - (n - 1)p_{\mathrm{in}}\bigr)\, t(1, 1) = 1 + (n - 1) p_{\mathrm{out}}\, t(2, 0). \tag{2}$$
   When $n = 1$, $n - 1 = 0$ and Eq. (2) directly yields $t(1, 1) = 1 = 1 \cdot H_1 + \frac{1 \cdot r (1 - 1)}{2}$.
   When $n \ge 2$, substituting Eq. (1) at $k = 2$, namely $t(2, 0) = \frac{1}{2 p_{\mathrm{in}}} + t(1, 1)$, into Eq. (2) gives
   $$\bigl(1 - (n - 1)(p_{\mathrm{in}} + p_{\mathrm{out}})\bigr)\, t(1, 1) = 1 + \frac{(n - 1) p_{\mathrm{out}}}{2 p_{\mathrm{in}}}.$$
   Using $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$, the left-hand coefficient is $1 - \frac{n-1}{n} = \frac{1}{n}$, so for all $n \ge 1$,
   $$t(1, 1) = n + \frac{n(n - 1)}{2}\, r = n H_1 + \frac{n r (n \cdot 1 - 1)}{1 + 1}. \tag{3}$$

3. For $2 \le k \le n$:
   $$\bigl(1 - (n - k)p_{\mathrm{in}}\bigr)\, t(k, 1) = 1 + k p_{\mathrm{in}}\, t(k - 1, 1) + k p_{\mathrm{out}}\, t(k, 0) + (n - k) p_{\mathrm{out}}\, t(k + 1, 0), \tag{4}$$
   where at $k = n$ the coefficient $(n - k) p_{\mathrm{out}}$ is zero, so $t(n+1, 0)$ does not enter. Substitute Eq. (1) for $t(k, 0)$ and (when $k < n$) for $t(k + 1, 0) = \frac{1}{(k+1)p_{\mathrm{in}}} + t(k, 1)$ into Eq. (4):
   $$\bigl(1 - (n - k)(p_{\mathrm{in}} + p_{\mathrm{out}})\bigr)\, t(k, 1) = 1 + \frac{p_{\mathrm{out}}}{p_{\mathrm{in}}} + \frac{n - k}{k + 1}\frac{p_{\mathrm{out}}}{p_{\mathrm{in}}} + k(p_{\mathrm{in}} + p_{\mathrm{out}})\, t(k - 1, 1).$$
   Using $p_{\mathrm{in}} + p_{\mathrm{out}} = \frac{1}{n}$ and $r = \frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}$, the coefficient of $t(k, 1)$ is $1 - \frac{n-k}{n} = \frac{k}{n}$, the coefficient of $t(k - 1, 1)$ is $\frac{k}{n}$, and the constant term on the right is
   $$1 + r\left(1 + \frac{n - k}{k + 1}\right) = 1 + \frac{n + 1}{k + 1}\, r.$$
   Multiplying both sides by $\frac{n}{k}$ gives the exact telescoping difference for every $2 \le k \le n$:
   $$t(k, 1) - t(k - 1, 1) = \frac{n}{k} + \frac{n(n + 1)r}{k(k + 1)} = \frac{n}{k} + n(n + 1)r\left(\frac{1}{k} - \frac{1}{k + 1}\right). \tag{5}$$

Summing Eq. (5) from $m = 2$ to $k$ ($1 \le k \le n$) and adding Eq. (3) yields
$$\begin{aligned}
t(k, 1) &= n + \frac{n(n - 1)}{2} r + n(H_k - 1) + n(n + 1)r\left(\frac{1}{2} - \frac{1}{k + 1}\right) \\
&= n H_k + n r \left(\frac{n - 1}{2} + \frac{n + 1}{2} - \frac{n + 1}{k + 1}\right) \\
&= n H_k + n r \left(n - \frac{n + 1}{k + 1}\right) \\
&= n H_k + \frac{n r (n k - 1)}{k + 1}.
\end{aligned}$$
Substituting this formula at $k - 1$ into Eq. (1), using $\frac{1}{p_{\mathrm{in}}} = n(1 + r)$ (since $n p_{\mathrm{in}}(1 + r) = 1$), gives for $2 \le k \le n$:
$$t(k, 0) = \frac{n(1 + r)}{k} + n H_{k-1} + \frac{n r (n(k - 1) - 1)}{k} = n H_{k-1} + \frac{n + n^2 r (k - 1)}{k}.$$
Finally, setting $k = n$ in $t(k, 1)$ gives $\frac{n r (n^2 - 1)}{n + 1} = n(n - 1) r$, and by Lemma 1,
$$T_N(p_{\mathrm{in}}, p_{\mathrm{out}}) = t(n, 1) = n H_n + n(n - 1)\frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}. \qquad \square$$

---

## 5. Corollaries: `bd`, `db`, and $\lambda$-mixed updating

**Corollary 1 (Proof of Brewster et al. 2025, Eq. (9): `bd` updating).** *Under neutral birth-death (`bd`) updating on the star graph $S_n$ with $N = n + 1$ vertices ($n \ge 1$) and no self-loops, starting from $N$ pairwise distinct types,*
$$T_N^{\mathrm{bd}} = n^3 - n^2 + n H_n.$$

*Proof.* Under `bd`, $p_{\mathrm{in}} = \frac{1}{Nn}$ and $p_{\mathrm{out}} = \frac{1}{N}$, so $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$ and $r = p_{\mathrm{out}}/p_{\mathrm{in}} = n$. By Theorem 1,
$$T_N^{\mathrm{bd}} = n H_n + n(n - 1)\cdot n = n^3 - n^2 + n H_n. \qquad \square$$

**Corollary 2 (Exact `db` consensus time on the star graph).** *Under neutral death-birth (`db`) updating on $S_n$ ($N = n + 1 \ge 2$, no self-loops), starting from $N$ pairwise distinct types,*
$$T_N^{\mathrm{db}} = n H_n + n - 1 = (N - 1) H_{N-1} + N - 2.$$

*Proof.* Under `db`, $p_{\mathrm{in}} = \frac{1}{N}$ and $p_{\mathrm{out}} = \frac{1}{Nn}$, so $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$ and $r = 1/n$. By Theorem 1,
$$T_N^{\mathrm{db}} = n H_n + n(n - 1)\cdot \frac{1}{n} = n H_n + n - 1. \qquad \square$$

**Corollary 3 (Exact $\lambda$-mixed consensus time on the star graph).** *Under neutral $\lambda$-mixed updating ($\lambda \in [0, 1]$, where each step is `bd` with probability $\lambda$ and `db` with probability $1-\lambda$) on $S_n$ ($N = n + 1 \ge 2$, no self-loops),*
$$T_N^{(\lambda)} = n H_n + n(n - 1)\,\frac{n\lambda + 1 - \lambda}{\lambda + n(1 - \lambda)}.$$
*In particular, for the unbiased half-half mixture $\lambda = 1/2$, $p_{\mathrm{in}} = p_{\mathrm{out}} = \frac{n+1}{2Nn} = \frac{1}{2n}$ ($r = 1$), giving*
$$T_N^{(1/2)} = n^2 - n + n H_n.$$

*Proof.* Under $\lambda$-mixed updating, $p_{\mathrm{in}}^{(\lambda)} = \frac{\lambda + n(1 - \lambda)}{Nn} > 0$ and $p_{\mathrm{out}}^{(\lambda)} = \frac{n\lambda + (1 - \lambda)}{Nn} > 0$, with $n(p_{\mathrm{in}}^{(\lambda)} + p_{\mathrm{out}}^{(\lambda)}) = 1$ and $r = \frac{n\lambda + 1 - \lambda}{\lambda + n(1 - \lambda)}$. Applying Theorem 1 completes the proof. $\square$

---

## 6. Structural remarks

1. **Identity $t(n, 0) = t(n, 1)$ under `bd` ($n \ge 2$).** From Theorem 1 at $k = n$ and $r = n$,
   $$t_{\mathrm{bd}}(n, 0) = n H_{n-1} + \frac{n + n^3(n - 1)}{n} = n H_{n-1} + 1 + n^3 - n^2 = n^3 - n^2 + n H_n = t_{\mathrm{bd}}(n, 1).$$
   Thus for $n \ge 2$, the expected time for the $n$ leaves alone to trace back to a single ancestor equals the expected time for all $n + 1$ vertices to trace back to a single ancestor. By contrast, at $n = 1$, $t_{\mathrm{bd}}(1, 0) = 0 \neq 1 = t_{\mathrm{bd}}(1, 1)$ because $(1, 0)$ is already a singleton.
2. **Comparison with the `db` coupon-collector bound.** Theorem 10 of Brewster et al. (2025, SI Section 9.4) notes the general `db` lower bound via the coupon-collector problem, requiring at least $N - 1 = n$ distinct death sites, which has expected time $n H_n$. Corollary 2 shows that on $S_n$ the exact `db` consensus time exceeds $n H_n$ by only $n - 1$ steps.
3. **Sensitivity to self-loops.** The conjecture and proof strictly assume a simple graph without self-loops. Adding a self-loop at the center (or at all vertices) changes $p_{\mathrm{in}}$ and $p_{\mathrm{out}}$ so that $n(p_{\mathrm{in}} + p_{\mathrm{out}}) < 1$, strictly increasing $T_N$ (confirmed explicitly in `checks/negative_controls.py` at $n = 3$).
