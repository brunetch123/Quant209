' Install Quant 209 for the current user. Double-click this file; nothing else is required.
'
' Installs to %LOCALAPPDATA%\Programs\Quant209, so no administrator rights are needed and
' nothing outside the user's own profile is touched. Creates a Desktop shortcut and a Start
' Menu entry. Run "Uninstall Quant 209.vbs" to remove all three.
'
' VBScript rather than a .bat file because a batch file flashes a console window, and the
' point of this installer is that a user never sees a terminal.

Option Explicit

Dim fso, shell, here, source, targetRoot, target, answer, desktopLink, menuLink

Set fso = CreateObject("Scripting.FileSystemObject")
Set shell = CreateObject("WScript.Shell")

here = fso.GetParentFolderName(WScript.ScriptFullName)
source = fso.BuildPath(here, "app")

If Not fso.FolderExists(source) Then
    MsgBox "Could not find the 'app' folder next to this installer." & vbCrLf & vbCrLf & _
           "Extract the whole ZIP first, then run this file from the extracted folder.", _
           vbCritical, "Quant 209 Setup"
    WScript.Quit 1
End If

targetRoot = fso.BuildPath(shell.ExpandEnvironmentStrings("%LOCALAPPDATA%"), "Programs")
target = fso.BuildPath(targetRoot, "Quant209")

If fso.FolderExists(target) Then
    answer = MsgBox("Quant 209 is already installed." & vbCrLf & vbCrLf & _
                    "Replace it with this version?", vbQuestion + vbYesNo, "Quant 209 Setup")
    If answer <> vbYes Then WScript.Quit 0
    On Error Resume Next
    fso.DeleteFolder target, True
    If Err.Number <> 0 Then
        MsgBox "Could not replace the existing installation." & vbCrLf & vbCrLf & _
               "Close Quant 209 if it is running, then try again.", vbCritical, "Quant 209 Setup"
        WScript.Quit 1
    End If
    On Error Goto 0
End If

If Not fso.FolderExists(targetRoot) Then fso.CreateFolder targetRoot

On Error Resume Next
fso.CopyFolder source, target, True
If Err.Number <> 0 Then
    MsgBox "Could not copy the program files." & vbCrLf & vbCrLf & Err.Description, _
           vbCritical, "Quant 209 Setup"
    WScript.Quit 1
End If
On Error Goto 0

' Desktop shortcut
Set desktopLink = shell.CreateShortcut(fso.BuildPath(shell.SpecialFolders("Desktop"), "Quant 209.lnk"))
desktopLink.TargetPath = fso.BuildPath(target, "Quant209.exe")
desktopLink.WorkingDirectory = target
desktopLink.Description = "PCB congener quantification"
desktopLink.IconLocation = fso.BuildPath(target, "Quant209.exe") & ", 0"
desktopLink.Save

' Start Menu entry
Set menuLink = shell.CreateShortcut(fso.BuildPath(shell.SpecialFolders("Programs"), "Quant 209.lnk"))
menuLink.TargetPath = fso.BuildPath(target, "Quant209.exe")
menuLink.WorkingDirectory = target
menuLink.Description = "PCB congener quantification"
menuLink.IconLocation = fso.BuildPath(target, "Quant209.exe") & ", 0"
menuLink.Save

answer = MsgBox("Quant 209 is installed." & vbCrLf & vbCrLf & _
                "A shortcut is on your Desktop and in the Start Menu." & vbCrLf & vbCrLf & _
                "Start it now?", vbInformation + vbYesNo, "Quant 209 Setup")

If answer = vbYes Then shell.Run """" & fso.BuildPath(target, "Quant209.exe") & """", 1, False
