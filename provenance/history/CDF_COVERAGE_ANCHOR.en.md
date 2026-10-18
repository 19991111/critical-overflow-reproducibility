# Selection-independent route for the four terminal-D CDF bounds

This proof distinguishes the CDF block from the post-hoc mean-betting block. It does not claim a new confirmatory micro replication or joint coverage for the whole follow-up workflow. Its numerical-coupling prerequisites are now supported by the received arithmetic supplement v2 under the stated arithmetic and ideal-word contracts; the current environment evidence is not a historical per-operation trace.

## 1. Reuse a genuinely earlier simultaneous event

The original v4 plan was fixed before its microscopic batch. For q=512 it declared eight individual empirical CDFs: two quantities (terminal D and terminal spill), lambda ∈ {1,2}, and the fixed/product kernels. It chose individual DKW width0.057 at1024 independent replicate labels and checked

16 exp[-2*1024*(0.057)^2] <0.025.

Let E_old be that simultaneous event. Then P(E_old^c)<=0.025 under the ideal model of iid uniform words across replicates and within every underlying stream. In particular, on E_old, all four terminal-D empirical CDFs differ from their population CDFs by at most0.057, and thus also by the looser0.060 used in the follow-up. This is an existing coverage event, not a post-hoc allocation to a newly chosen stake. It holds at every x simultaneously. Reusing consequences of the same event does not require paying for it again each time a consequence is stated.

## 2. Condition only the new reference on the old history

Let H be the entire old-data/design history available before the new reference draws. Conditional on H, the two reference targets, coefficients, sample size32768, numerical scheme, error bounds and seeds/ideal stream contract were fixed. Under independent new ideal reference replicate streams, the two empirical reference CDFs obey width0.012 bands with conditional failure probability bounded by

4 exp[-2*32768*(0.012)^2] <0.0004.

The last scalar bound follows from the already checked rational Taylor lower bound exp[2*32768*(0.012)^2]>10000. The reference bound may depend on the old design history through its fixed conditional target, but the probability bound is uniform in that history. Therefore averaging conditional probabilities gives P(E_ref^c)<=0.0004. No independence between E_old and E_ref is required for the union bound.

Thus P(E_old intersect E_ref)>=1-0.025-0.0004=0.9746. This is a CDF-only statement. It is not the error budget for all old experiments, all later studies or the post-hoc mean procedure.

## 3. Map the event to the unchanged interval calculations

For each lambda ∈ {1,2}, the deterministic population-level SDE coupling gives

F_comp(x-r)-tau <= F_SDE(x) <= F_comp(x+r)+tau,

where (r,tau)=(0.010,0.005) or(0.020,0.005). Tau is a population CDF approximation allowance, not a confidence-event failure probability. The Euler, coefficient, driver and arithmetic assumptions must justify this population inequality. On E_ref, replace F_comp with its empirical CDF and add0.012. On E_old, replace the micro population CDF with its empirical CDF and add at most0.060. Taking the exact directed suprema over x gives precisely the earlier recorded intervals, with their0.060+0.012+0.005 expansion. No endpoints or thresholds are recomputed here.

Because these inequalities hold simultaneously over all x and the same four predeclared marginals, subsequent inspection of x, selection of a subset of these four comparisons, or comparison of the resulting intervals with the stated tolerance does not invalidate this event. The scientific study remains a post-hoc diagnostic. This reasoning does not extend the event to arbitrary new marginals, scales, transformations or mean-betting recipes.

## 4. What remains unproved

The original micro mean betting intervals were designed after seeing the old data and lack the same pre-data coverage anchor. The CDF event above does not imply uniformity of product capitals over retrospectively chosen stakes or clipping rules. Consequently this CDF-only route cannot justify the former claim that the entire eight-comparison selection workflow has joint coverage at least95%. The received v2 arithmetic supplement supports the numerical prerequisites under its contract. Current environment checks support that contract, while the original dynamic register state was not logged and actual arithmetic remains a conditional premise.
