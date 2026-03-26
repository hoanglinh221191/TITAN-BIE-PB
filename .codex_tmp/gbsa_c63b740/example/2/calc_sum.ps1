$file = $args[0]
if (-not $file) { Write-Error "No file provided"; exit 1 }
$content = Get-Content $file | Select-Object -Skip 1
$sum = 0
foreach ($line in $content) {
    if ([string]::IsNullOrWhiteSpace($line)) { continue }
    $parts = $line.Trim() -split '\s+'
    if ($parts.Length -ge 3) {
        try {
            $v1 = [double]$parts[1]
            $v2 = [double]$parts[2]
            $sum += $v1 + $v2
        } catch {}
    }
}
Write-Host "Sum for $file : $sum"
