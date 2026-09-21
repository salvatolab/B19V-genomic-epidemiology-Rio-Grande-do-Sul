"""Testes estatisticos dos dados epidemiologicos da Tabela 1 (Gregianini 2026).

Entrada: contagens agregadas da Tabela 1 do manuscrito (nao ha microdado aqui).
Saida: os valores que vao para a nova secao de analise estatistica.
"""
import numpy as np
from scipy import stats

YEARS = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024]
TESTED = [206, 564, 522, 764, 973, 432, 164, 311, 54, 128]
POS    = [58, 109, 128, 73, 110, 38, 6, 16, 0, 15]
FEM    = [37, 70, 93, 43, 60, 22, 3, 8, 0, 9]
MAL    = [21, 39, 35, 30, 50, 16, 3, 8, 0, 6]
AGE = {"0-6m": 26, "7m-1y": 84, "2-5y": 68, "6-10y": 62, "11-20y": 70,
       "21-30y": 103, "31-40y": 92, "41-50y": 29, "51-60y": 17, "60y+": 2}

def wilson(k, n):
    lo, hi = stats.binomtest(k, n).proportion_ci(method="wilson")
    return 100 * lo, 100 * hi

n_tested, n_pos = sum(TESTED), sum(POS)
n_fem, n_mal = sum(FEM), sum(MAL)
print(f"testadas {n_tested}  positivas {n_pos}  ({100*n_pos/n_tested:.1f}%, IC95 {wilson(n_pos, n_tested)[0]:.1f}-{wilson(n_pos, n_tested)[1]:.1f})")

# 1. proporcao de mulheres entre os positivos
bt = stats.binomtest(n_fem, n_fem + n_mal, 0.5)
print(f"\nsexo feminino {n_fem}/{n_fem+n_mal} = {100*n_fem/(n_fem+n_mal):.1f}% "
      f"(IC95 {wilson(n_fem, n_fem+n_mal)[0]:.1f}-{wilson(n_fem, n_fem+n_mal)[1]:.1f})  "
      f"binomial bilateral p = {bt.pvalue:.3g}")

# 2. positividade por ano
tab = np.array([POS, [t - p for t, p in zip(TESTED, POS)]])
chi2, p, dof, exp = stats.chi2_contingency(tab)
print(f"\npositividade por ano: chi2 = {chi2:.1f}, gl = {dof}, p = {p:.3g}  "
      f"(menor esperado {exp.min():.1f})")
for y, t, po in zip(YEARS, TESTED, POS):
    print(f"   {y}: {po}/{t} = {100*po/t:.1f}%  (IC95 {wilson(po,t)[0]:.1f}-{wilson(po,t)[1]:.1f})")

# 3. razao de sexos ao longo dos anos (2023 sem casos fica fora)
idx = [i for i in range(len(YEARS)) if FEM[i] + MAL[i] > 0]
tab2 = np.array([[FEM[i] for i in idx], [MAL[i] for i in idx]])
chi2b, pb, dofb, expb = stats.chi2_contingency(tab2)
print(f"\nsexo por ano ({len(idx)} anos): chi2 = {chi2b:.1f}, gl = {dofb}, p = {pb:.3g}  "
      f"(menor esperado {expb.min():.1f})")
simul = stats.chi2_contingency(tab2, correction=False,
                              method=stats.MonteCarloMethod(rng=42, n_resamples=20000))
print(f"   mesmo teste por permutacao (esperados pequenos): p = {simul.pvalue:.3g}")

# 4. faixas etarias
tot = sum(AGE.values())
k = AGE["21-30y"] + AGE["31-40y"]
print(f"\n21-40 anos: {k}/{tot} = {100*k/tot:.1f}% (IC95 {wilson(k,tot)[0]:.1f}-{wilson(k,tot)[1]:.1f})")
obs = np.array(list(AGE.values()))
chi2c, pc = stats.chisquare(obs)
print(f"distribuicao etaria uniforme entre as 10 faixas: chi2 = {chi2c:.1f}, gl = {len(obs)-1}, p = {pc:.3g}")
