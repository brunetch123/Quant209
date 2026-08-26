' Remove Quant 209 for the current user. Double-click this file.
'
' Deletes the program folder from %LOCALAPPDATA%\Programs\Quant209 and both shortcuts.
' Nothing else was ever written, so nothing else needs cleaning up -- in particular the
' analysis data and any reports you exported are left exactly where they are.

Option Explicit

Dim fso, shell, target, answer, desktopLink, menuLink, failed

Set fso = CreateObject("Scripting.FileSystemObject")
Set shell = CreateObject("WScript.Shell")

target = fso.BuildPath(shell.ExpandEnvironmentStrings("%LOCALAPPDATA%"), "Programs\Quant209")
desktopLink = fso.BuildPath(shell.SpecialFolders("Desktop"), "Quant 209.lnk")
menuLink = fso.BuildPath(shell.SpecialFolders("Programs"), "Quant 209.lnk")

If Not fso.FolderExists(target) And Not fso.FileExists(desktopLink) Then
    MsgBox "Quant 209 does not appear to be installed for this user.", _
           vbInformation, "Quant 209 Setup"
    WScript.Quit 0
End If

answer = MsgBox("Remove Quant 209?" & vbCrLf & vbCrLf & _
                "Your data and any exported reports are not affected.", _
                vbQuestion + vbYesNo, "Quant 209 Setup")
If answer <> vbYes Then WScript.Quit 0

failed = False
On Error Resume Next
If fso.FolderExists(target) Then
    fso.DeleteFolder target, True
    If Err.Number <> 0 Then failed = True
    Err.Clear
End If
If fso.FileExists(desktopLink) Then fso.DeleteFile desktopLink, True : Err.Clear
If fso.FileExists(menuLink) Then fso.DeleteFile menuLink, True : Err.Clear
On Error Goto 0

If failed Then
    MsgBox "Some files could not be removed." & vbCrLf & vbCrLf & _
           "Close Quant 209 if it is running, then run this again.", vbExclamation, "Quant 209 Setup"
Else
    MsgBox "Quant 209 has been removed.", vbInformation, "Quant 209 Setup"
End If
