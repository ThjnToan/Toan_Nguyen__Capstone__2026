$f = "C:\Users\Laptop K1\OneDrive\Documents\MATLAB\main.tex"
$t = [System.IO.File]::ReadAllText($f, [System.Text.Encoding]::UTF8)

# Add Z row to calibration table
$old = [string]::Concat("			Battery Price (SS) & ", '$', "\bar{P}_{bat}", '$', " & 1.0 & Normalized \\")
$new = [string]::Concat("			Battery Price (SS) & ", '$', "\bar{P}_{bat}", '$', " & 1.0 & Normalized \\`r`n			TFP Scale & ", '$', "Z", '$', " & 1.131 & Calibrated so ", '$', "Y_{ss}=1", '$', " (sole output normalization); absorbs residual scaling from ", '$', "L_{ss}=\tfrac{1}{3}", '$', ", ", '$', "\bar{E}=0.15", '$', " \\")
$t = $t.Replace($old, $new)

# Update steady state table note (just the most distinctive part)
$t = $t.Replace("73.2\% reflects the small open economy", "73.4\% reflects the small open economy")
$t = $t.Replace("`$\tau_{ss} = 1.65\%`$ is endogenously determined from the proportional tax rule `$\tau_t = I_{grid,t}/Y_t`$. Output `$Y = 0.861`$ is the endogenous level emerging from CES aggregation of value-added and the energy input `$\bar{E} = 0.15`$; its numerical value reflects the chosen unit system (labor time endowment normalized so `$L_{ss} = \tfrac{1}{3}`$) and carries no additional normalization.", "`$\tau_{ss}=1.43\%`$ equals `$\delta_g K_{g,ss}/Y_{ss}`$; with `$Y_{ss}=1`$ this simplifies to `$\tau_{ss}=\delta_g K_{g,ss}`$ exactly. Output `$Y_{ss}=1`$ is enforced by the TFP scale `$Z=1.131`$ (sole normalization).")

[System.IO.File]::WriteAllText($f, $t, [System.Text.Encoding]::UTF8)
Write-Host "Done"
if ($t -match "TFP Scale") { Write-Host "Z row: YES" } else { Write-Host "Z row: NO" }
