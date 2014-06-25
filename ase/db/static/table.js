function open_row(row, id, tid) {
    var table = document.getElementById('rows');
    r = row.rowIndex;
    if (row.classList.contains('open')) {
        row.classList.remove('open');
        table.deleteRow(r + 1);
    } else {
        request = new XMLHttpRequest();
        request.open('GET', '/open_row/' + id + '?x=' + tid, true);
        request.onload = function() {
            data = request.responseText;
            var newrow = table.insertRow(r + 1);
            row.classList.add('open');
            var cell = newrow.insertCell(0);
            cell.colSpan = 100;
            cell.innerHTML = data;
        }
        request.send();
    }
}
