$env:PATH += ';C:\dynare\6.5\preprocessor'
$pinfo = New-Object System.Diagnostics.ProcessStartInfo
$pinfo.FileName = 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe'
$pinfo.Arguments = '-batch "cd(''C:/Users/Laptop K1/OneDrive/Documents/MATLAB''); run_dynare_proptax"'
$pinfo.RedirectStandardOutput = $true
$pinfo.RedirectStandardError = $true
$pinfo.UseShellExecute = $false
$p = New-Object System.Diagnostics.Process
$p.StartInfo = $pinfo
$p.Start() | Out-Null
$stdout = $p.StandardOutput.ReadToEnd()
$stderr = $p.StandardError.ReadToEnd()
$p.WaitForExit()
($stdout + $stderr) | Out-File 'C:/Users/Laptop K1/OneDrive/Documents/MATLAB/matlab_output.txt' -Encoding utf8
exit $p.ExitCode
