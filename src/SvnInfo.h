#ifndef SVN_INFO_INCLUDED
#define SVN_INFO_INCLUDED

const char* VersionInfo::WCREV = "$WCREV$";//Replaced with the highest commit revision in the working copy
const char* VersionInfo::WCDATE = "$WCDATE$";//Replaced with the commit date/time of the highest commit revision. By default, international format is used: yyyy-mm-dd hh:mm:ss. Alternatively, you can specify a custom format which will be used with strftime(), for example: $WCDATE=%a %b %d %I:%M:%S %p$
const char* VersionInfo::WCNOW = "$WCNOW$";//Replaced with the current system date/time. This can be used to indicate the build time. Time formatting is as described above. 
const char* VersionInfo::WCRANGE = "$WCRANGE$";//Replaced with the update revision range in the working copy. If the working copy is in a consistent state, this will be a single revision. If the working copy contains mixed revisions, either due to being out of date, or due to a deliberate update-to-revision, then the range will be shown in the form 100:200 
const char* VersionInfo::WCMIXED = "$WCMIXED?Mixed update revision:Not mixed$";//$WCMIXED?TText:FText$ is replaced with TText if there are mixed update revisions, or FText if not. 
const char* VersionInfo::WCMODS = "$WCMODS?Modified:Not modified$";//$WCMODS?TText:FText$ is replaced with TText if there are local modifications, or FText if not. 
const char* VersionInfo::WCURL = "$WCURL$";//Replaced with the repository URL of the working copy path passed to SubWCRev.

#endif