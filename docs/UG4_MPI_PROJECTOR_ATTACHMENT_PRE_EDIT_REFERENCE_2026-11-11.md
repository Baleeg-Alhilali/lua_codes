# UG4 MPI projector-attachment pre-edit reference — 2026-11-11

This note records the exact recoverable state immediately before testing the
suspected MPI projector-attachment defect. It is a rollback reference, not a
claim that the local modification is correct.

## Repository identities

- UG4 core repository: `/Users/alhilaba/UG4_promesh_ogrid/ugcore`
- UG4 core `HEAD`: `46b9b57a8c8398939bf1b6a4661396c863f6a723`
- NavierStokes repository: `/Users/alhilaba/UG4_promesh_ogrid/plugins/NavierStokes`
- NavierStokes `HEAD`: `8df0ce228312397f10638219ba39383dbf59b4ac`
- `navier_stokes_fv1.cpp` SHA-256 before this test:
  `a423b74aa06874b01f33a419b563a2ddeb82e55c0b0a2d9f8fa002639200a12b`

The NavierStokes FV1 source itself was clean relative to its repository at this
checkpoint. The relevant local change is in UG4 core's
`ugbase/lib_grid/parallelization/distribution.cpp`.

## Exact relevant pre-edit diff from UG4-core HEAD

```diff
@@ -1734,9 +1734,14 @@ void AddGlobalAttachmentsToSerializer (
 		GridDataSerializationHandler& handler,
 		Grid& grid)
 {
-	const vector<string>& attachmentNames = GlobalAttachments::declared_attachment_names();
+	// Declaration order may differ between ranks (e.g. when a serialized
+	// projector declares its attachments while being reconstructed).  Every
+	// rank must nevertheless install serializers in exactly the same order.
+	vector<string> attachmentNames = GlobalAttachments::declared_attachment_names();
+	std::sort(attachmentNames.begin(), attachmentNames.end());
 	for(size_t i = 0; i < attachmentNames.size(); ++i){
 		const string& name = attachmentNames[i];
+		if(name == "npSurfParams") continue; // diagnostic: isolate custom attachment transfer
 		if(GlobalAttachments::is_attached<TElem>(grid, name)){
 			GlobalAttachments::add_data_serializer<TElem>(handler, grid, name);
 		}
@@ -1752,7 +1757,10 @@ static void SynchronizeAttachedGlobalAttachments (
 	UG_DLOG(LG_DIST, 1, "SynchronizeAttachedGlobalAttachments start\n");
 
 //	now make sure that the same global attachments are attached everywhere
-	const vector<string>&  names = GlobalAttachments::declared_attachment_names();
+	// Compare and attach by a canonical name order.  Equal counts alone are not
+	// sufficient because later serializer order is communication-significant.
+	vector<string> names = GlobalAttachments::declared_attachment_names();
+	std::sort(names.begin(), names.end());
```

To restore the immediate pre-test state after the proposed one-line edit,
reinsert this line immediately after `const string& name = attachmentNames[i];`:

```cpp
if(name == "npSurfParams") continue; // diagnostic: isolate custom attachment transfer
```

Do not reset or check out the whole UG4-core repository: it contains other
unrelated local work that must be preserved.
